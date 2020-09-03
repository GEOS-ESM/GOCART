#!/usr/bin/env python

"""
  Python wrapper for ana_lpe.x.

  Chem_Ext3d.py aer_f aod_d 

"""

import os
import sys

from optparse  import OptionParser   # Command-line args  

from gfio import GFIO
from MAPL import config

import warnings
warnings.filterwarnings("ignore")

#---------------------------------------------------------------------

if __name__ == "__main__":

    expid = "a0000"
    dir   = '.'

#   System Dependent defaults
#   -------------------------
    if os.path.exists("/nobackup/1/ARCTAS"):  # calculon
        dir = '/nobackup/1/ARCTAS/'
        aer = dir + 'Y%y4/M%m2/d5_arctas_02.inst3d_aer_v.%y4%m2%d2_%h2%n2z.nc'
        ext = "Y%y4/M%m2/%s.ext.sfc.%y4%m2%d2_%h2%n2z.nc"
    elif os.path.exists('/discover/nobackup/projects/gmao/iesa/'): # Discover
        dir = '/discover/nobackup/projects/gmao/iesa/aerosol/data/ARCTAS/'
        aer = dir + 'Y%y4/M%m2/D%d2/d5_arctas_02.inst3d_aer_v.%y4%m2%d2_%h2%n2z.nc4'
        ext = "Y%y4/M%m2/%s.ext.sfc.%y4%m2%d2_%h2%n2z.nc"
    else: # Generic default
        aer =  "%s.aer_Nv.eta.%y4%m2%d2_%h2%n2z.nc"
        ext =  "%s.ext.sfc.%y4%m2%d2_%h2%n2z.nc"
        

    rc_file = "ext.rc"
    Nx = "2"
    Ny = "4"
    channel = "532"

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [options] nymd nhms",
                          version='1.0.0' )

    parser.add_option("-c", "--channel", dest="channel", default=channel,
                      help="Wavelength in nm (default=%s)"%channel )

    parser.add_option("-i", "--aer", dest="aer", default=aer,
                      help="input  aerosol mixing ratio file template (default=%s)"%aer )

    parser.add_option("-o", "--ext", dest="ext", default=ext,
                      help="output extinction file template (default=%s)"%ext )

    parser.add_option("-r", "--rcfile", dest="rc", default=rc_file,
                      help="resource file (default=%s)"%rc_file )
    
    parser.add_option("-X", "--Nx", dest="Nx", default=Nx,
                      help="number of PEs to decompose longitude (default=%s)"%Nx )

    parser.add_option("-Y", "--Ny", dest="Ny", default=Ny,
                      help="number of PEs to decompose latitude (default=%s)"%Ny )
    
    parser.add_option("-x", "--expid", dest="expid", default=expid,
                      help="experiment Id (default=%s)"%expid )
    
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose")

    options, args = parser.parse_args()
    
    if len(args) < 2:
        parser.error("not enough input arguments")
    else:
        nymd = args[0]
        nhms = args[1]

#   Expand file name templates
#   --------------------------
    options.aer = config.strTemplate(options.aer,expid=options.expid,
                                       nymd=nymd,nhms=nhms)
    options.ext = config.strTemplate(options.ext,expid=options.expid,
                                       nymd=nymd,nhms=nhms)

#   Get file dimensions
#   -------------------
    f = GFIO(options.aer,'r')

#   Load rc file
#   ------------
    cf = config.Config(options.rc)

#   Update rc file with user specified parameters
#   ---------------------------------------------
    cf('ExpId',options.expid)
    cf('Layout_Nx',options.Nx)
    cf('Layout_Ny',options.Ny)
    cf('IM_World',f.im)
    cf('JM_World',f.jm)
    cf('LM_World',f.km)
    cf('nymd',nymd)
    cf('nhms',nhms)
    cf('BANDS',"%sE-9"%options.channel)
    if options.verbose:
        cf('verbose',".TRUE.")
    else:
        cf('verbose',".FALSE.")
    cf('aer_filename',options.aer)
    cf('ext_filename',options.ext)
    cf.save(rcfile=options.rc) # save updated rc file

    nPE = int(options.Nx) * int(options.Ny) 

#   Run the Fortran binary
#   ----------------------
    rc = os.system("mpirun -np %d ext_calculator.xx"%nPE)
    if rc:
        raise RuntimeError, "rc=%d on return from 'ext_calculator.xx'"%rc
