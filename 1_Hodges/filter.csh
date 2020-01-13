#!/bin/csh

#SBATCH --account=ccl
#SBATCH --job-name=TrackingT42
#SBATCH -c 1
#SBATCH --time=00:30:00
#SBATCH -a 0-39:1

 set year = 1979
 @ year=$year + $SLURM_ARRAY_TASK_ID
 set HOME = ~/TRACK-ERAI
 set WORK = `pwd`

 set target = /rigel/ccl/users/yy2865/Output/ERA-Interim/N128/850/T42
 
 cd $HOME/bin

cat > rsplice.$year.opt << EOZ
n
0
4
n
1
y
n
g
n
1
170
1
71
y
1
n
0
1
${target}/${year}/objout.new.A.$year
${target}/${year}/tdump.A.$year
1
y
y
3
1
8
84
n
n
n
n
a
y
n
n
0
EOZ
   nice ./track.A.$year < rsplice.$year.opt
   mv -f $HOME/outdat/*A.${year}* ${target}/${year}
   
   mkdir ${targer}/${year}/script
   mv $WORK/*.csh ${target}/${year}/script
   mv $WORK/*.in  ${target}/${year}/script
   mv $WORK/ymp.${year} ${target}/${year}/script
  
# End of script
