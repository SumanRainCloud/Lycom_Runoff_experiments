#!/bin/bash

################################
# Settings
################################
#
# name of binary
#
vname="25Jan_100species_30year_veg_current" #t3fO  #_21"  #"H6C6sp"
expname="25Jan_100species_30year_veg_current" #-n8c36
#
model="lycom_${vname}.x"
#
# initial number of species
#
nspecies=100
#
# simulation length
#
simlength=30  # years
startout=1
endout=30
#
# output interval
#
intvout=24 #168 #24     # hours
#
# switch (0/1) for submission to queuing system
#
lsf=1
#
# switches
#
lBSC=".false."
lNOHONO=".false."
lnoVeg=".true."
#
# PATHS
#
if [[ $lsf == 1 ]]
then
  basedir=/beegfs/u/bas7785/lycom_current

  bindir=${basedir}/binaries

  rundir=/beegfs/u/bas7785/simulations/${expname}

  inputdir=/beegfs/u/bas7785/postprocessed/
 
else
  basedir=/beegfs/u/bas7785/lycom_current

  bindir=${basedir}/binariesL

  rundir=${basedir}/simulations/${expname}
  #rundir=/home/philipp/work/supervision/elevatedCO2_bryophytes/simulations/${expname}

  inputdir=${basedir}/
fi
#
# lsf settings
#
ptn="std"  	# [ compute / compute2 ]
ntasks=24        # nodes ???
cpus_per_task=8          # max, tasks p. node [ (24 / 36) ] ???
#
runID=${expname} #${rundir##*/}
#
################################

# make run directory

if [ -e ${rundir} ]
then
  echo "run directory already exists!"
  exit
fi

mkdir ${rundir}

# copy/link input files to run directory

if [[ $lsf == 1 ]]
then

  ln -s ${inputdir}/srad_new.nc ${rundir}/srad.nc
  ln -s ${inputdir}/lrad_new.nc ${rundir}/lrad.nc
  ln -s ${inputdir}/tair_new.nc ${rundir}/tair.nc
  ln -s ${inputdir}/rhum_new.nc ${rundir}/rhum.nc
  ln -s ${inputdir}/rain_new.nc ${rundir}/rain.nc
  ln -s ${inputdir}/snow_new.nc ${rundir}/snow.nc
  ln -s ${inputdir}/wind_new.nc ${rundir}/wind.nc
  ln -s ${inputdir}/lycom_specpar ${rundir}/lycom_specpar
  ln -s ${inputdir}/landsea_new.nc ${rundir}/landsea.nc
  ln -s ${inputdir}/biome_new.nc ${rundir}/biome.nc
  ln -s ${inputdir}/ETrmax.nc ${rundir}/ETrmax.nc
  ln -s ${inputdir}/ETdavg.nc ${rundir}/ETdavg.nc
  ln -s ${inputdir}/LAI_new.nc ${rundir}/LAI.nc
  ln -s ${inputdir}/SAI_new.nc ${rundir}/SAI.nc
  ln -s ${inputdir}/SSA_new.nc ${rundir}/SSA.nc
  if [[ \${lNOHONO} == .true. ]]; then ln -s ${inputdir}/NO_HONO_fSatBins ${rundir}/NO_HONO_fSatBins ; fi
else

  cd ${inputdir}

  cp srad lrad tair rhum rain snow wind libry_specpar LAI SAI SSA biome ${rundir}

  if [[ $lNOHONO == .true. ]]
  then
    cp NO_HONO_fSatBins ${rundir}
  fi
fi

cp ${bindir}/binary_${vname}/${model} ${rundir}

# start simulation

cd ${rundir}

cat > lycom_namelist << EOF
&lycompar
cyear0=2010,
tsindata0=1,
lastyear=${simlength},
tsl=3600,
tpos0=1,
yearout1=${startout},
yearoutX=${endout},
outint=${intvout},
nSites=1,
p_nspec=${nspecies},
frac_s_init=0.0,
specout=.true.,
interCan=.false.,
NOHONO=${lNOHONO},
BSCtypes=${lBSC},
noVeg=${lnoVeg},
inDirect=0
&end
EOF
#outdiurnal=.true.,

if [[ $lsf == 1 ]]      # parallel run
then

cat > lycom_batch << EOF
#!/bin/bash
#SBATCH --job-name=run_${runID}         # Specify job name
#SBATCH --partition=${ptn}              # Specify partition name (compute/compute2)
#SBATCH --nodes=1                   # Specify number of nodes
#SBATCH --ntasks=24        # Specify number of tasks per node (max. 24/36)
#SBATCH --cpus-per-task=8
#SBATCH --time="12:00:00"
#SBATCH --mail-type=FAIL                # Notify user by email in case of job failure
#SBATCH --mail-user=suman.halder@uni-hamburg.de	# email address
#SBATCH --export=NONE					# ?
#SBATCH --output=job_${runID}.o%j       # File name for standard output
#SBATCH --error=job_${runID}.e%j        # File name for standard error output

source /sw/batch/init.sh

# Environment settings

#module switch env env/2015Q2-intel16-impi
module switch env env/2023Q4-gcc-openmpi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# launch program
echo "launching program" ${model}
#mpirun ./${model}

mpirun --bind-to core ./${model}
EOF

  # submit script to cluster

  #sbatch --partition=${ptn} lycom_batch
  sbatch lycom_batch
else

  ./${model}


  cp ${basedir}/run_lycom.sh ${rundir}
fi


exit

# rm srad.nc
# rm lrad.nc
# rm tair.nc
# rm rhum.nc
# rm rain.nc
# rm snow.nc
# rm wind.nc
# rm libry_specpar
# rm landsea.nc
# rm biome.nc
# rm LAI.nc
# rm SAI.nc
# rm SSA.nc
# if [[ \${lNOHONO} == .true. ]]; then rm NO_HONO_fSatBins; fi

