#!/bin/csh -f
#
# Computes aerosol extention files. See usage: below for details.
#

echo "Beginning Exctinction Calculator  $argv"

# set up executable names (assumes they are on the path)
  set chemaod = "${FVROOT}/bin/Chem_Aod3d.x"
  set lcv2prs = "${FVROOT}/bin/lcv2prs.x"

# Did I pass in a filename to operate on?
  if ($#argv < 1) then
   goto usage
  endif

# Check options
  if ( "$1" == "-clean" ) then
    set clean = 1
    shift
  else
    set clean = 0
  endif

# Setup the input filename
  set inpfile = $1

# Parse the name of the input filename
  set expid      = `echo $inpfile:r:r:r`
  set inpfiletyp = `echo $inpfile:r:r:e`
  set datetag    = `echo $inpfile:r:e`

# Trickery to pull out the YYYYMMDD and HHMMSS from filename
  set datevalid = `echo $datetag   | cut -d"+" -f2`
  set YYYYMMDD  = `echo $datevalid | cut -c1-8`
  set HHMMSS    = `echo $datevalid | cut -c10-13`00

# Check the input file passed
# If it is not of type "filetyp" then exit
  if ( "$inpfiletyp" =~ *3d_aer_v* ) then
   echo "Starting extinction calculation for file $inpfile"
  else
   exit 0 # Not an error, exit gracefully
  endif
  /bin/rm -f tau3d.nc4
# Now run the AOD calculator
# With no other options will produce output file tau3d.nc4
   set cmd = "$chemaod -t Aod_CALIPSO.rc $inpfile"     
   echo $cmd
   $cmd
  if ( $status ) then
    echo $0": error running the extinction calculator"
    exit 2
  endif

# rename the output file
  set extvfile = `echo $inpfile | sed -e 's/_aer_v/_ext_v/g'`
  \mv -f tau3d.nc4 $extvfile

# Now do the pressure level interpolation
  set extpfile = `echo $inpfile | sed -e 's/_aer_v/_ext_p/g'`
  set savfile = `echo $inpfile | sed -e 's/_aer_v/_sav_v/g'`
  $lcv2prs -nStep=1 -date $YYYYMMDD -start $HHMMSS \
           -vars @tavg3d_ext_p -rc tavg3d_ext_p.rc \
           -o $extpfile $extvfile
  if ( $status ) then
    echo $0": error on return from lcv2prs - aborting"
    exit 3
  endif

# Optional clean-up: delete input and intermediate files
  if ( $clean ) then

#        /bin/cp $inpfile $extvfile $extpfile /explore/nobackup/dao_ops/colarco
         /bin/mv $inpfile $savfile
        \rm -f $inpfile $extvfile 
  endif

  echo $0":All done"
# Kludge for MPI run
 if ( $?MPI_ENVIRONMENT ) ${FVROOT}/bin/makeiau.x
 exit

  exit 0

#--------------------------------------------------------------
usage:

cat <<EOF
NAME
     extinction_calculator - Calculates 3D extinction from aer_v files
          
SYNOPSIS
     extinction_calculator.csh [-clean] aer_v_filename
           
DESCRIPTION
     This script will run the extinction calculator Chem_Aod3d.x
     against the input file, generating an interim 3d extinction file
     named EXPID.tavg3d_ext_v.DATETAG.nc4.  The interim 3d extinction
     file is then run through lcv2prs.x to turn into a pressure level
     file.  The input file is assumed to conform to the following
     dile anme convention:
                EXPID.FTYPE.DATETAG.nc4
     where

     EXPID   is the experiment id (e.g., d5_arctas_02)
     FTYPE   the calculator will accept the following file types:
                inst3d_aer_v 
                tavg3d_aer_v
                inst3d_aer_vunz
                tavg3d_aer_vunz
     DATETAG is the forecast start and valid time date string in the
             format YYYYMMDD_HHz+YYYYMMDD_HHMMz

    The final pressure level file is named exactly as the input
    file with the string "aer_v" replaced by "ext_p".

    When option "-clean" is specified, both the input file as well
    as the intermediate "ext_v" file will be deleted.

REQUIRES
    It requires the the following executables to be on your path:
       Chem_Aod3d.x
       lcv2prs.x
    The following resource files should be on your current directory:
       Aod_CALIPSO.rc
       tavg3d_ext_p.rc

BUGS
    It does not work with assimilation files (will need to revise the
    DATETAG handling).

AUTHOR
    Written by Peter.R.Colarco@nasa.gov 
    Revised by Arlindo.daSilva@nasa.gov

EOF

exit 1
