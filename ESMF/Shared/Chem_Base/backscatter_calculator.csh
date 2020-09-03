#!/bin/csh -f
#
# Computes aerosol backscatter files files. See usage: below for details.
#

# set up executable names (assumes they are on the path)
  set bin = Chem_Aod3d.x

# Did I pass in a filename to operate on?
  if ($#argv < 2) then
   goto usage
  endif

# Setup the input filename
  set aerFile = "$1"
  shift

  if ( ! (-e $aerFile) ) then
    echo $0": cannot find input aer_v file $aerFile"
    exit 1
  endif

  foreach ch ( $argv )

     set chnm = ${ch}nm
     set rcFile = Aod3d_$chnm.rc
     set absFile = ( `echo $aerFile | sed -e "s/_aer_v/_abs-${chnm}_v/g"` )
  
     if ( ! (-e $rcFile) ) then
        echo $0": cannot find input rc file $rcFile"
        exit 1
     endif

     /bin/rm -f $absFile

     set cmd = "$bin -t $rcFile -o $absFile $aerFile"

     echo $cmd
     $cmd

  end

  echo $0": All done"
  exit 0

#--------------------------------------------------------------
usage:

cat <<EOF
NAME
     backscatter_calculator - Calculates 3D backscatter from aer_v files
          
SYNOPSIS
     backscatter_calculator.csh aer_v_filename channel(s)
           
DESCRIPTION
     This script will run the 3D aerosol "calculator" Chem_Aod3d.x
     against the input file, generating a 3D output file for each
     wavelength. On input, "channel" refers to the desired wavenumbers, 
     usually one or more among (355nm,532nm,1064nm).

EXAMPLE

    % backscatter_calculator.csh \
          d5_arctas_02.inst3d_aer_v.20080629_1200z.hdf 355 532
                                
    This will create 2 files:

       d5_arctas_02.inst3d_abs-355nm_v.20080629_1200z.hdf 
       d5_arctas_02.inst3d_abs-532nm_v.20080629_1200z.hdf 

    The following resource files are required:

       Chem_MieRegistry.rc
       Abs_355nm.rc
       Abs_532nm.rc

AUTHOR
    Arlindo.daSilva@nasa.gov

EOF

exit 1
