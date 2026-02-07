#!/bin/bash

# SETTINGS
#-------------------------------

para=1

#-------------------------------



a0=25Jan_100species_30year_veg_current

if [ -e "lycom.x" ]
then
  rm lycom.x
fi

if [ -e "lycom*.mod" ]
then
  rm lycom*.mod
fi

if [[ $para == 1 ]]
then

  # DKRZ
  #source /sw/rhel6-x64/etc/profile.mistral 	 	# make module command available
  #module load intel/18.0.2 intelmpi/2018.1.163
  #module add intel/18.0.2 intelmpi/2018.1.163
  
  #mpiifort -fp-model strict -free -real-size 64 -O0 \
  #         $(/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/bin/nf-config --fflags --flibs) \
  #         -g -traceback -ftrapuv -debug all -gen-interfaces \
  #         libry_par.f90 libry_common.f90 libry_global.f90 libry_nolocal.f90 libry_noopt.f90 libry.f90 libry_main.f90 -o libry.x
  #         #libry_par.f90 libry_common.f90 libry_global.f90 libry_local.f90 libry_opt.f90 libry.f90 libry_main.f90 -o libry.x

  ########

  # Hummel-2
  source /sw/modules/init.sh

  module switch env env/2023Q4-gcc-openmpi
  
  mpifort -ffree-form -fdefault-real-8 -ffree-line-length-none \
      $(nf-config --fflags --flibs) \
    -fbacktrace -fcheck=all -fbounds-check -fimplicit-none \
    -Og -g -ffpe-trap=invalid,zero,overflow,underflow \
     lycom_par.f90 lycom_common.f90 lycom_global.f90 lycom_nolocal.f90 \
     lycom_noopt.f90 lycom.f90 lycom_land.f90 lycom_interface.f90 lycom_main.f90 -o lycom.x

           #libry_par.f90 libry_common.f90 libry_global.f90 libry_local.f90 libry_opt.f90 libry.f90 libry_main.f90 -o libry.x

  # -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interfaces -diag-enable=sc \
else

  gfortran -ffree-line-length-256 \
           lycom_par.f90 lycom_common.f90 lycom_noglobal.f90 lycom_local.f90 lycom_noopt.f90 lycom.f90 lycom_land.f90 lycom_interface.f90 lycom_main.f90 -o lycom.x

  # -fcheck=all 
fi

rm lycom*.mod

if [[ $para == 1 ]]
then
  mkdir -p binaries  
  if [ -e "binaries/binary_${a0}" ]
  then
    if [[ $1 == "c" ]]
    then
      rm -r ./binaries/binary_${a0}
   
    else
      echo "binary directory exists!"
      exit
    fi
  fi

  mkdir binaries/binary_${a0}

  mv lycom.x binaries/binary_${a0}/lycom_${a0}.x
  cp lycom*f90 binaries/binary_${a0}
  cp compile_lycom.sh binaries/binary_${a0}

else
  mkdir -p binariesL
  if [ -e "binariesL/binary_${a0}" ]
  then
    if [[ $1 == "c" ]]
    then
      rm -r binariesL/binary_${a0}
    else
      echo "binary directory exists!"
      exit
    fi
  fi

  mkdir binariesL/binary_${a0}

  mv lycom.x binariesL/binary_${a0}/lycom_${a0}.x
  cp lycom*f90 binariesL/binary_${a0}
  cp compile_lycom.sh binariesL/binary_${a0}

fi

exit
