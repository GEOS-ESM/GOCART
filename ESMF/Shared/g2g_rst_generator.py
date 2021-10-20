# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Written by E. Sherman on 05/05/2021

# DESCRIPTION: This program will generate the necessary GOCART2G restart files
#              if provided a legacy GOCART restart file (gocart_internal_rst).
#              Every use-case has likely not been tested. If you have any issues
#              please contact Elliot Sherman at elliot.m.sherman@nasa.gov.


# HOW TO USE:  This script is meant to be run on NCCS Discover. Please load the
#              following module "python/GEOSpyD/Ana2019.03_py3.7" to run.
#              After loading the module, type the following in your terminal:
#              "python3 g2g_rst_generator.py" and follow the prompts.
#              
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from netCDF4 import Dataset

# prompt user for gocart restart file
# -----------------------------------
gocart_rst = input("Path to gocart_internal_rst: ")

# ask user if brown carbon is included in restart
# ------------------------------------------------
answer = input("Is brown carbon included in the gocart_internal_rst? (yes/no): ")

if answer.upper() == 'YES':
  aeroList = ['du','ss','ni','su','caoc','cabc','cabr']
if answer.upper() == 'NO':
  aeroList = ['du','ss','ni','su','caoc','cabc']

# read original gocart restart
# ----------------------------
dataOrig = Dataset(gocart_rst, mode='r')

# function to create GOCART2G restarts
# ------------------------------------
def createAeroRestart (dataOrig, aeroName):

  # create new restart file
  # -----------------------
  data = Dataset(aeroName + '_internal_rst', mode='w', format='NETCDF4')

  # create dimensions
  # -----------------
  if aeroName == 'du' or aeroName == 'ss':
    unknown_dim1 = data.createDimension("unknown_dim1", None)
  lat = data.createDimension("lat", dataOrig.dimensions['lat'].size)
  lev = data.createDimension("lev", dataOrig.dimensions['lev'].size)
  lon = data.createDimension("lon", dataOrig.dimensions['lon'].size)
  time = data.createDimension("time", dataOrig.dimensions['time'].size)

  # create new lon variable
  # -----------------------
  lon = data.createVariable("lon", "d", ("lon"))
  lon.long_name = dataOrig.variables['lon'].long_name
  lon.units = dataOrig.variables['lon'].units
  # fill new lon with original lon's values
  lon = data.variables['lon']
  lon[:] = dataOrig.variables['lon'][:]

  # create new lat variable
  # -----------------------
  lat = data.createVariable("lat", "d", ("lat"))
  lat.long_name = dataOrig.variables['lat'].long_name
  lat.units = dataOrig.variables['lat'].units
  # fill new lat with old lat's values
  lat = data.variables['lat']
  lat[:] = dataOrig.variables['lat'][:]

  # create new lev variable
  # ----------------------
  lev = data.createVariable("lev", "d", ("lev"))
  lev.coordinate = dataOrig.variables['lev'].coordinate
  lev.formulaTerms = dataOrig.variables['lev'].formulaTerms
  lev.long_name = dataOrig.variables['lev'].long_name
  lev.positive = dataOrig.variables['lev'].positive
  lev.standard_name = dataOrig.variables['lev'].standard_name
  lev.units = dataOrig.variables['lev'].units
  # fill new lev with original lev's values
  lev = data.variables['lev']
  lev[:] = dataOrig.variables['lev'][:]

  # create new time variable
  # -------------------------
  time = data.createVariable("time", "d", ("time"))
  time.units = dataOrig.variables['time'].units
  time.long_name = "time"
  # fill new time with original time's values
  time = data.variables['time']
  time[:] = dataOrig.variables['time'][:]


  # create new aerosol variable and fill with original aerosol values
  # -----------------------------------------------------------------
  if aeroName == 'du':
    aero1 = data.createVariable("DU", "f4", ("unknown_dim1","lev", "lat", "lon"))
    aero1.long_name = "Dust Mixing Ratio All bins"
    aero1.units = "kg kg-1"

    aero1 = data.variables['DU']
    aero1[0,:,:,:] = dataOrig.variables['du001'][:]
    aero1[1,:,:,:] = dataOrig.variables['du002'][:]
    aero1[2,:,:,:] = dataOrig.variables['du003'][:]
    aero1[3,:,:,:] = dataOrig.variables['du004'][:]
    aero1[4,:,:,:] = dataOrig.variables['du005'][:]

  elif aeroName == 'ss':
    aero1 = data.createVariable("SS", "f4", ("unknown_dim1","lev", "lat", "lon"))
    aero1.long_name = "Sea Salt Mixing Ratio All bins"
    aero1.units = "kg kg-1"

    aero1 = data.variables['SS']
    aero1[0,:,:,:] = dataOrig.variables['ss001'][:]
    aero1[1,:,:,:] = dataOrig.variables['ss002'][:]
    aero1[2,:,:,:] = dataOrig.variables['ss003'][:]
    aero1[3,:,:,:] = dataOrig.variables['ss004'][:]
    aero1[4,:,:,:] = dataOrig.variables['ss005'][:]

  elif aeroName == 'ni':
    aero1 = data.createVariable("NH3", "f4", ("lev", "lat", "lon"))
    aero1.long_name = "Ammonia (NH3, gas phase)"
    aero1.units = "kg kg-1"
    aero2 = data.createVariable("NH4a", "f4", ("lev", "lat", "lon"))
    aero2.long_name = "Ammonium ion (NH4+, aerosol phase)"
    aero2.units = "kg kg-1"
    aero3 = data.createVariable("NO3an1", "f4", ("lev", "lat", "lon"))
    aero3.long_name = "Nitrate size bin 001"
    aero3.units = "kg kg-1"
    aero4 = data.createVariable("NO3an2", "f4", ("lev", "lat", "lon"))
    aero4.long_name = "Nitrate size bin 002"
    aero4.units = "kg kg-1"
    aero5 = data.createVariable("NO3an3", "f4", ("lev", "lat", "lon"))
    aero5.long_name = "Nitrate size bin 003"
    aero5.units = "kg kg-1"

    aero1 = data.variables['NH3']
    aero1[:,:,:] = dataOrig.variables['NH3'][:]
    aero2 = data.variables['NH4a']
    aero2[:,:,:] = dataOrig.variables['NH4a'][:]
    aero3 = data.variables['NO3an1']
    aero3[:,:,:] = dataOrig.variables['NO3an1'][:]
    aero4 = data.variables['NO3an2']
    aero4[:,:,:] = dataOrig.variables['NO3an2'][:]
    aero5 = data.variables['NO3an3']
    aero5[:,:,:] = dataOrig.variables['NO3an3'][:]

  elif aeroName == 'su':
    aero1 = data.createVariable("DMS", "f4", ("lev", "lat", "lon"))
    aero1.long_name = "Dimethylsulphide"
    aero1.units = "kg kg-1"
    aero2 = data.createVariable("MSA", "f4", ("lev", "lat", "lon"))
    aero2.long_name = "Methanesulphonic acid"
    aero2.units = "kg kg-1"
    aero3 = data.createVariable("SO2", "f4", ("lev", "lat", "lon"))
    aero3.long_name = "Sulphur dioxide"
    aero3.units = "kg kg-1"
    aero4 = data.createVariable("SO4", "f4", ("lev", "lat", "lon"))
    aero4.long_name = "Sulphate aerosol"
    aero4.units = "kg kg-1"

    aero1 = data.variables['DMS']
    aero1[:,:,:] = dataOrig.variables['DMS'][:]
    aero2 = data.variables['MSA']
    aero2[:,:,:] = dataOrig.variables['MSA'][:]
    aero3 = data.variables['SO2']
    aero3[:,:,:] = dataOrig.variables['SO2'][:]
    aero4 = data.variables['SO4']
    aero4[:,:,:] = dataOrig.variables['SO4'][:]

  elif aeroName == 'caoc':
    aero1 = data.createVariable("CAphilicCA.oc", "f4", ("lev", "lat", "lon"))
    aero1.long_name = "Hydrophilic Carbonaceous Aerosol"
    aero1.units = "kg kg-1"
    aero2 = data.createVariable("CAphobicCA.oc", "f4", ("lev", "lat", "lon"))
    aero2.long_name = "Hydrophobic Carbonaceous Aerosol"
    aero2.units = "kg kg-1"

    aero1 = data.variables['CAphilicCA.oc']
    aero1[:,:,:] = dataOrig.variables['OCphilic'][:]
    aero2 = data.variables['CAphobicCA.oc']
    aero2[:,:,:] = dataOrig.variables['OCphobic'][:]

  elif aeroName == 'cabc':
    aero1 = data.createVariable("CAphilicCA.bc", "f4", ("lev", "lat", "lon"))
    aero1.long_name = "Hydrophilic Carbonaceous Aerosol"
    aero1.units = "kg kg-1"
    aero2 = data.createVariable("CAphobicCA.bc", "f4", ("lev", "lat", "lon"))
    aero2.long_name = "Hydrophobic Carbonaceous Aerosol"
    aero2.units = "kg kg-1"

    aero1 = data.variables['CAphilicCA.bc']
    aero1[:,:,:] = dataOrig.variables['BCphilic'][:]
    aero2 = data.variables['CAphobicCA.bc']
    aero2[:,:,:] = dataOrig.variables['BCphobic'][:]

  elif aeroName == 'cabr':
    aero1 = data.createVariable("CAphilicCA.br", "f4", ("lev", "lat", "lon"))
    aero1.long_name = "Hydrophilic Carbonaceous Aerosol"
    aero1.units = "kg kg-1"
    aero2 = data.createVariable("CAphobicCA.br", "f4", ("lev", "lat", "lon"))
    aero2.long_name = "Hydrophobic Carbonaceous Aerosol"
    aero2.units = "kg kg-1"

    aero1 = data.variables['CAphilicCA.br']
    aero1[:,:,:] = dataOrig.variables['BRCphilic'][:]
    aero2 = data.variables['CAphobicCA.br']
    aero2[:,:,:] = dataOrig.variables['BRCphobic'][:]

  # Close file
  data.close()

# Run the function and generate a restart for each aerosol
# --------------------------------------------------------
for aeroName in aeroList:
  createAeroRestart(dataOrig, aeroName)



