
#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=20:00:00
#PBS -j eo 
#PBS -q opterons
#PBS -N PatientCNT6Case0
cd /mnt/io1/compbio/home/aghos/workdir/prostate_testosterone_scan/SimDirCNT6/Case0
echo "Job started on `hostname` at `date`"
echo " "
python combined_ode_boolean.py -p Model_Input.json &> Output.log
echo "Job Ended on `hostname` at `date`"
echo " "
