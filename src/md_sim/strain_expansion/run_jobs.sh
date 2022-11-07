#/bin/sh
yourfilenames=`ls ./lammps.*`
STEP=0

for slurmfile in $yourfilenames
do
   temp=${slurmfile//[^0-9]/}
   echo $temp
   squeue -u bh326
   echo $slurmfile
   sbatch $slurmfile
   mv $slurmfile done/

   STEP=$((STEP+1))
   if [ "$STEP" -eq "8" ]
      then
      echo "Found three files"
      squeue -u bh326
      break
   fi
done
#sbatch $slurmfile
#mv $eachfile $slurmfile inputs/done/
