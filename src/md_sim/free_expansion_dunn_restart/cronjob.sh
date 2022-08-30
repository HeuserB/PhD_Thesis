#!/bin/sh
ssh -tt Titan /bin/bash <<'EOT'
squeue -u bh326 > log
MAXJOBS=2
RUNNINGJOBS=$(wc --lines < log)
let RUNNINGJOBS-=1
rm log
echo "Currently running $RUNNINGJOBS of $MAXJOBS jobs."

cd free_expansion_dunn/

if [ $RUNNINGJOBS < $MAXJOBS ]; then
    COUNT=$(ls inputs_free_expansion/lammps.* | wc -l)
    if [ COUNT > 0 ]; then
        echo "Free space and $COUNT lammps files found"
        sh free_expansion_dunn/run_more.sh
    else
        echo "No more lammps files found. Testing if all jobs have finished."
        LOGDIR=/data/bh326/LAMMPS/free_expansion_dunn
        python inputs_free_expansion_restart/file_comparison.py $LOGDIR
    fi
fi

exit
EOT