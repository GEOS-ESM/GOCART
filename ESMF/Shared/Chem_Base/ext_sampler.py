#!/usr/bin/env python

"""
    Utility to compute optical properties along a sampled file.

    Ed Nowottnick, January, 2015.

"""

import os
import MieObs_
from types import *
from netCDF4 import Dataset
from mieobs import VNAMES, getAOPext, getAOPscalar
from numpy import zeros, arange, array, ones, zeros, interp, isnan, ma, NaN, squeeze, transpose, shape, asarray
from math import pi, sin, cos, asin, acos

from types           import *
from optparse        import OptionParser
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser
from csv             import DictReader

from MAPL           import Config, eta
from MAPL.constants import *
from pyobs          import ICARTT
from pyobs.sgp4     import getTrack 

# Generic Lists of Varnames and Units
VNAMES_DU = ['DU001','DU002','DU003','DU004','DU005']
VNAMES_SS = ['SS001','SS002','SS003','SS004','SS005']
VNAMES_BC = ['BCPHOBIC','BCPHILIC']
VNAMES_OC = ['OCPHOBIC','OCPHILIC']
VNAMES_SU = ['SO4']
MieVarsNames = ['ext','scatext','backscat','aback_sfc','aback_toa','depol','ext2back','tau','ssa','g']
MieVarsUnits = ['km-1','km-1','km-1 sr-1','sr-1','sr-1','unitless','sr','unitless','unitless','unitless']

class MieCalc(object):
    pass
    """                                                                                  
    Generic container for Variables
    """

#---
def getVars(inFile):
    """
#    Parse input file, create variable dictionary
#    """

    Vars       = MieCalc()
    file       = Dataset(inFile)
    names      = file.variables.keys()
    MIENAMES   = names
    for n, name in enumerate(MIENAMES):
        var = file.variables[name]
        if (name == 'trjLon') or (name == 'stnLon'):
            name = 'LONGITUDE'
        if (name == 'trjLat') or (name == 'stnLat'):
            name = 'LATITUDE'
        name = name.upper()
        size = len(var.shape)
        if size == 3:
            setattr(Vars,name,var[:,:,:])
        if size == 2:
            setattr(Vars,name,var[:,:])
        if size == 1:
            setattr(Vars,name,var[:])
    return Vars        

#---
def computeMie(Vars, channel, varnames, rcFile):
    """
#    Computes optical quantities and combines into a dictionary
#   """
    
    #STN Sampled?
    if options.station:
        NAMES = VNAMES + ['PS','DELP','RH','AIRDENS']
        nstn = len(Vars.STATION)
        nobs = len(Vars.TIME)
        for v in range(nstn):
            VarsIn = MieCalc()
            for n, name in enumerate(NAMES):
                Var = getattr(Vars,name)
                size = len(Var.shape)
                if (size == 2):
                    #1D Variables, ex. PS
                    setattr(VarsIn,name,Var[v,:])
                if size == 3:
                    #2D Variables, ex. DU001
                    setattr(VarsIn,name,Var[v,:,:])
  
            if (v==0):
                tau,ssa,g = getAOPscalar(VarsIn,channel,vnames=varnames,vtypes=varnames,Verbose=True,rcfile=rcFile)
                ext,sca,backscat,aback_sfc,aback_toa,depol = getAOPext(VarsIn,channel,I=None,vnames=varnames,vtypes=varnames,Verbose=True,rcfile=rcFile)
                ext2back = ext/backscat
                MieVars = {"ext":[ext],"scatext":[sca],"backscat":[backscat],"aback_sfc":[aback_sfc],"aback_toa":[aback_toa],"depol":[depol],"ext2back":[ext2back],"tau":[tau],"ssa":[ssa],"g":[g]}
            else:
                tau,ssa,g = getAOPscalar(VarsIn,channel,vnames=varnames,vtypes=varnames,Verbose=True,rcfile=rcFile)
                ext,sca,backscat,aback_sfc,aback_toa,depol = getAOPext(VarsIn,channel,I=None,vnames=varnames,vtypes=varnames,Verbose=True,rcfile=rcFile)
                ext2back = ext/backscat
                MieVars['ext'].append(ext)
                MieVars['scatext'].append(sca)
                MieVars['backscat'].append(backscat)
                MieVars['aback_sfc'].append(aback_sfc)
                MieVars['aback_toa'].append(aback_toa)
                MieVars['depol'].append(depol) 
                MieVars['ext2back'].append(ext2back)
                MieVars['tau'].append(tau)
                MieVars['ssa'].append(ssa)
                MieVars['g'].append(g)

    #TRJ Sampled?
    else:
        tau,ssa,g = getAOPscalar(Vars,channel,vnames=varnames,vtypes=varnames,Verbose=True,rcfile=rcFile)
        ext,sca,backscat,aback_sfc,aback_toa,depol = getAOPext(Vars,channel,I=None,vnames=varnames,vtypes=varnames,Verbose=True,rcfile=rcFile)
        ext2back = ext/backscat
        MieVars = {"ext":[ext],"scatext":[sca],"backscat":[backscat],"aback_sfc":[aback_sfc],"aback_toa":[aback_toa],"depol":[depol],"ext2back":[ext2back],"tau":[tau],"ssa":[ssa],"g":[g]}

    return MieVars

#---
def writeNC ( stations, lons, lats, tyme, isotimeIn, MieVars, MieVarsNames,
              MieVarsUnits, inFile, outFile, zlib=False):
    """
    Write a NetCDF file with sampled GEOS-5 variables along the satellite track
    described by (lon,lat,tyme).
    """
    km = 72

    # Open NC file
    # ------------
    nc = Dataset(outFile,'w',format=options.format)

    # Set global attributes
    # ---------------------
    nc.title = 'GEOS-5 Sampled Aerosol Optical Properties File'
    nc.institution = 'NASA/Goddard Space Flight Center'
    nc.source = 'Global Model and Assimilation Office'
    nc.history = 'Created from sampled GEOS-5 collections'
    nc.references = 'n/a'
    nc.comment = 'This file contains sampled GEOS-5 aerosol optical properties.'
    nc.contact = 'Ed Nowottnick <edward.p.nowottnick@nasa.gov>'
    nc.Conventions = 'CF'
    nc.inFile = inFile
 
    # Create dimensions
    # -----------------
    nt = nc.createDimension('time',len(tyme))
    if options.station:
        ns = nc.createDimension('station',len(stations))
    ls = nc.createDimension('ls',19)
    if km>0:
        nz = nc.createDimension('lev',km)
    x = nc.createDimension('x',1)
    y = nc.createDimension('y',1)

    if options.station:
        # Station names
        # -------------
        stnName_ = nc.createVariable('stnName','S1',('station','ls'),zlib=zlib)
        stnName_.long_name = 'Station Names'
        stnName_.axis = 'e'
        stnName_[:] = stations[:]   

    # Coordinate variables
    # --------------------
    time = nc.createVariable('time','i4',('time',),zlib=zlib)
    time.long_name = 'Time'
    t0 = tyme[0]
    isot0 = isotimeIn[0]
    date0 = ''.join(isot0[:10])
    time0 = ''.join(isot0[-8:])
    time.units = 'seconds since '+date0+' '+time0
    time[:] = tyme
    if km > 0: # pressure level not supported yet
        lev = nc.createVariable('lev','f4',('lev',),zlib=zlib)
        lev.long_name = 'Vertical Level'
        lev.units = 'km'
        lev.positive = 'down'
        lev.axis = 'z'
        lev[:] = range(1,km+1)

    # Add fake dimensions for GrADS compatibility
    # -------------------------------------------
    x = nc.createVariable('x','f4',('x',),zlib=zlib)
    x.long_name = 'Fake Longitude for GrADS Compatibility'
    x.units = 'degrees_east'
    x[:] = zeros(1)
    y = nc.createVariable('y','f4',('y',),zlib=zlib)
    y.long_name = 'Fake Latitude for GrADS Compatibility'
    y.units = 'degrees_north'
    y[:] = zeros(1)
    if options.station:
        e = nc.createVariable('station','i4',('station',),zlib=zlib)
        e.long_name = 'Station Ensemble Dimension'
        e.axis = 'e'
        e.grads_dim = 'e'
        e[:] = range(len(stations))
    
    # Lat/Lon Coordinates
    # ----------------------
    if options.station:
        lon = nc.createVariable('longitude','f4',('station',),zlib=zlib)
        lon.long_name = 'Longitude'
        lon.units = 'degrees_east'
        lon[:] = lons[:]
        lat = nc.createVariable('latitude','f4',('station',),zlib=zlib)
        lat.long_name = 'Latitude'
        lat.units = 'degrees_north'
        lat[:] = lats[:]
    else:
        lon = nc.createVariable('longitude','f4',('time',),zlib=zlib)
        lon.long_name = 'Longitude'
        lon.units = 'degrees_east'
        lon[:] = lons[:]
        lat = nc.createVariable('latitude','f4',('time',),zlib=zlib)
        lat.long_name = 'Latitude'
        lat.units = 'degrees_north'
        lat[:] = lats[:]        
    
    # Time in ISO format if so desired
    # ---------------------------------
    isotime = nc.createVariable('isotime','S1',('time','ls'),zlib=zlib)
    isotime.long_name = 'Time (ISO Format)'
    isotime[:] = isotimeIn[:]
      
    # Write each variable
    # --------------------------------------------------
    for n, name in enumerate(MieVarsNames):

        var = squeeze(MieVars[name])
        size = len(var.shape)
        if options.station:
            if size == 2:
                var = asarray([var])
                size = len(var.shape)
        if size == 3:
            dim = ('station','time','lev')
        if size == 2:
            dim = ('time','lev')
        if size == 1:
            dim = ('time')
        this = nc.createVariable(name,'f4',dim,zlib=zlib)
        this.standard_name = name
        this.units = MieVarsUnits[n]
        this.missing_value = MAPL_UNDEF
        if options.station:
            this[:] = transpose(var,(0,2,1))
        else:
            this[:] = transpose(var)

    # Close the file
    # --------------
    nc.close()

    if options.verbose:
        print " <> wrote %s file %s"%(options.format,options.outFile)
    
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    
    format = 'NETCDF3_CLASSIC'
    inFile  = 'trj_sampler.nc'
    outFile = 'ext_sampler.nc'
    channel = (532,)
    rcFile = 'Aod3d_532nm.rc'

#   Parse command line options
#   --------------------------
    parser = OptionParser()

    parser.add_option("-i", "--input", dest="inFile", default=inFile,
              help="Sampled input file")

    parser.add_option("-o", "--output", dest="outFile", default=outFile,
              help="Output file containing optical properties")

    parser.add_option("-r", "--rc", dest="rcFile", default=rcFile,
              help="Resource file pointing to optical tables")

    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC, NETCDF3_64BIT or EXCEL (default=%s)"%format )

    parser.add_option("-c", "--channel", dest="channel", default=channel,
              help="Channel for Mie calculation")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode")

    parser.add_option("--du",
                      action="store_true", dest="dust",
                      help="Dust Only")

    parser.add_option("--ss",
                      action="store_true", dest="seasalt",
                      help="Seasalt Only")

    parser.add_option("--su",
                      action="store_true", dest="sulfate",
                      help="Sulfate Only")

    parser.add_option("--bc",
                      action="store_true", dest="bcarbon",
                      help="Black Carbon Only")

    parser.add_option("--oc",
                      action="store_true", dest="ocarbon",
                      help="Organic Carbon Only")

    parser.add_option("--stn",
                      action="store_true", dest="station",
                      help="Input File is from stn_sampler.py")

    (options, args) = parser.parse_args()

         
    # Create consistent file name extension
    # -------------------------------------
    name, ext = os.path.splitext(options.outFile)
    if ext.upper() == '.XLS':
        options.format = 'EXCEL'
    if 'NETCDF4' in options.format:
        options.outFile = name + '.nc4'
    elif 'NETCDF3' in options.format:
        options.outFile = name + '.nc'
    elif 'EXCEL' in options.format:
        options.outFile = name + '.xls'
    else:
        raise ValueError, 'invalid extension <%s>'%ext
 

    # Get Variables
    # --------------------------
    Vars = getVars(options.inFile)

    # Run Mie Calculator and Write Output Files
    # --------------------------
    if options.station:
        StnNames = Vars.STNNAME
    else:
        StnNames = ''

    channelIn = float(options.channel)
    MieVars = computeMie(Vars,channelIn,VNAMES,options.rcFile)
    writeNC(StnNames,Vars.LONGITUDE,Vars.LATITUDE,Vars.TIME,Vars.ISOTIME,MieVars,MieVarsNames,MieVarsUnits,options.inFile,options.outFile)

    if options.dust:
        outFile = options.outFile+'.dust'
        MieVars = computeMie(Vars,channelIn,VNAMES_DU,options.rcFile)
        writeNC(StnNames,Vars.LONGITUDE,Vars.LATITUDE,Vars.TIME,Vars.ISOTIME,MieVars,MieVarsNames,MieVarsUnits,options.inFile,outFile)

    if options.seasalt:
        outFile = options.outFile+'.ss'
        MieVars = computeMie(Vars,channelIn,VNAMES_SS,options.rcFile)
        writeNC(StnNames,Vars.LONGITUDE,Vars.LATITUDE,Vars.TIME,Vars.ISOTIME,MieVars,MieVarsNames,MieVarsUnits,options.inFile,outFile)

    if options.sulfate:
        outFile = options.outFile+'.su'
        MieVars = computeMie(Vars,channelIn,VNAMES_SU,options.rcFile)
        writeNC(StnNames,Vars.LONGITUDE,Vars.LATITUDE,Vars.TIME,Vars.ISOTIME,MieVars,MieVarsNames,MieVarsUnits,options.inFile,outFile)

    if options.bcarbon:
        outFile = options.outFile+'.bc'
        MieVars = computeMie(Vars,channelIn,VNAMES_BC,options.rcFile)
        writeNC(StnNames,Vars.LONGITUDE,Vars.LATITUDE,Vars.TIME,Vars.ISOTIME,MieVars,MieVarsNames,MieVarsUnits,options.inFile,outFile)

    if options.ocarbon:
        outFile = options.outFile+'.oc'
        MieVars = computeMie(Vars,channelIn,VNAMES_OC,options.rcFile)
        writeNC(StnNames,Vars.LONGITUDE,Vars.LATITUDE,Vars.TIME,Vars.ISOTIME,MieVars,MieVarsNames,MieVarsUnits,options.inFile,outFile)
   
