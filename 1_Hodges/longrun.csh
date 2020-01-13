#!/bin/csh

#SBATCH --account=ccl
#SBATCH --job-name=TrackingT42
#SBATCH -c 2
#SBATCH --time=3:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yujia@ldeo.columbia.edu
#SBATCH -a 0-39:1

     set year = 1979
     @ year=$year + $SLURM_ARRAY_TASK_ID
     set HOME = ~/TRACK-ERAI
     set WORK = `pwd`

     set target = /rigel/ccl/users/yy2865/Output/ERA-Interim/N128/850/T42
     set source = /rigel/ccl/users/yy2865/ERA-Interim/N128/850/AMJJASO/T42

     # create a directory for outputs
     mkdir -p ${target}/${year}

     # soft link the data to TRACK working directory
     cd $HOME/indat
     ln -s $source/$year.Asia.nc $year.Asia.nc
     
     cd $HOME/bin
     
     nice ./track.A.$year < $WORK/optlst.Hodges.in
     mv -f $HOME/outdat/*A.${year}* ${target}/${year}

     mkdir ${targer}/${year}/script
     mv $WORK/*.csh ${target}/${year}/script
     mv $WORK/*.in  ${target}/${year}/script
     mv $WORK/ymp.${year} ${target}/${year}/script
  
# End of script
