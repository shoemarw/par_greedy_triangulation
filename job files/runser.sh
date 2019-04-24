#!/bin/bash
#
#SBATCH --job-name=needOnly40mins
#SBATCH --output=ser_tri.txt
#SBATCH --ntasks=1
echo "ya-hoooo!"
for x in 707 1000 1414 2000 2828; do
	echo "== Serial version $x points =="
	srun "tri" "test"$x"pts93" | grep "Gent"
done
