#!/usr/bin/env python

"""
  Python wrapper for ana_lpe.x.

  ana_lde.py aer_f aod_d 

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
        dir_f = '/nobackup/1/ARCTAS/'
        dir_a = '/home/adasilva/GAAS/%s/chem/'
        aer_f = dir_f + 'Y%y4/M%m2/d5_arctas_02.inst3d_aer_v.%y4%m2%d2_%h2%n2z.nc'
        aer_a = dir_a + "Y%y4/M%m2/%s.aer_a.eta.%y4%m2%d2_%h2%n2z.nc4"
        aod_f = dir_a + "Y%y4/M%m2/%s.aod_f.sfc.%y4%m2%d2_%h2%n2z.nc"
        aod_d = dir_a + "Y%y4/M%m2/%s.aod_d.sfc.%y4%m2%d2_%h2%n2z.nc"
    elif os.path.exists('/discover/nobackup/projects/gmao/iesa/'): # Discover
        dir_f = '/discover/nobackup/projects/gmao/iesa/aerosol/data/ARCTAS/'
        dir_a = '/discover/nobackup/projects/gmao/iesa/aerosol/experiments/GAAS/%s/chem/'
        aer_f = dir_f + 'Y%y4/M%m2/D%d2/d5_arctas_02.inst3d_aer_v.%y4%m2%d2_%h2%n2z.nc4'
        aer_a = dir_a + "Y%y4/M%m2/D%d2/%s.aer_a.eta.%y4%m2%d2_%h2%n2z.nc4"
        aod_f = dir_a + "Y%y4/M%m2/%s.aod_f.sfc.%y4%m2%d2_%h2%n2z.nc"
        aod_d = dir_a + "Y%y4/M%m2/%s.aod_d.sfc.%y4%m2%d2_%h2%n2z.nc"
    else: # Generic default
        aer_f =  "%s.aer_f.eta.%y4%m2%d2_%h2%n2z.nc"
        aer_a =  "%s.aer_a.eta.%y4%m2%d2_%h2%n2z.nc"
        aod_f =  "%s.aod_f.sfc.%y4%m2%d2_%h2%n2z.nc"
        aod_d =  "%s.aod_d.sfc.%y4%m2%d2_%h2%n2z.nc"
        

    rc_file = "lde.rc"
    Nx = "2"
    Ny = "4"

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [options] nymd nhms",
                          version='1.0.0' )

    parser.add_option("-A", "--aer_a", dest="aer_a", default=aer_a,
                      help="output aerosol concentration analysis   file template (default=%s)"%aer_a )

    parser.add_option("-F", "--aer_f", dest="aer_f", default=aer_f,
                      help="input  aerosol concentration background file template (default=%s)"%aer_f )

    parser.add_option("-D", "--dir", dest="dir", default=dir,
                      help="directory name to append to file names (default=%s)"%dir )

    parser.add_option("-d", "--aod_d", dest="aod_d", default=aod_d,
                      help="output AOD analysis increment file template (default=%s)"%aod_d )

    parser.add_option("-f", "--aod_f", dest="aod_f", default=aod_f,
                      help="output AOD background file template (default=%s)"%aod_f )

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
    options.aer_f = config.strTemplate(options.aer_f,expid=options.expid,
                                       nymd=nymd,nhms=nhms)
    options.aer_a = config.strTemplate(options.aer_a,expid=options.expid,
                                       nymd=nymd,nhms=nhms)
    options.aod_f = config.strTemplate(options.aod_f,expid=options.expid,
                                       nymd=nymd,nhms=nhms)
    options.aod_d = config.strTemplate(options.aod_d,expid=options.expid,
                                       nymd=nymd,nhms=nhms)

#   Append directory
#   ----------------
    if options.aer_f[0] not in ('/','.'):
        options.aer_f = options.dir+'/'+options.aer_f
    if options.aer_a[0] not in ('/','.'):
        options.aer_a = options.dir+'/'+options.aer_a
    if options.aod_f[0] not in ('/','.'):
        options.aod_f = options.dir+'/'+options.aod_f
    if options.aod_d[0] not in ('/','.'):
        options.aod_d = options.dir+'/'+options.aod_d

#   Get file dimensions
#   -------------------
    f = GFIO(options.aer_f,'r')
    d = GFIO(options.aod_d,'r')

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
    cf('CM_World',d.km)
    cf('nymd',nymd)
    cf('nhms',nhms)
    if options.verbose:
        cf('verbose',".TRUE.")
    else:
        cf('verbose',".FALSE.")
    cf('aer_ana_filename',options.aer_a)
    cf('aer_bkg_filename',options.aer_f)
    cf('aod_bkg_filename',options.aod_f)
    cf('aod_inc_filename',options.aod_d)
    cf.save(rcfile="lde.rc") # save updated rc file

    nPE = int(options.Nx) * int(options.Ny) 

#   Run the Fortran binary
#   ----------------------
    rc = os.system("mpirun -np %d ana_lde.x"%nPE)
    if rc:
        raise RuntimeError, "rc=%d on return from 'ana_lde.x'"%rc
