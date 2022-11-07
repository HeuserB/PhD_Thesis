#!/bin/sh
squeue -u bh326 > log
OPENJOBS=$(sed -n '$=' log)
OPENJOBS=$((OPENJOBS-1))
OPENJOBS=$((8-OPENJOBS))

echo $OPENJOBS
rm log

yourfilenames=`ls ./lammps.*`
STEP=0

for slurmfile in $yourfilenames
do
   if [ "$STEP" -eq "$OPENJOBS" ]
      then
      printf "No More Open Jobs\nCurrently running:\n"
      squeue -u bh326
      break
   fi

   STEP=$((STEP+1))

   #squeue -u bh326
   echo $slurmfile
   sbatch $slurmfile
   mv $slurmfile done/
done
#sbatch $slurmfile
#mv $eachfile $slurmfile inputs/done/

