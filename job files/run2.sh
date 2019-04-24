#!/bin/bash
#
#SBATCH --job-name=icanhazclustertiem
#SBATCH --output=par_tri.txt
#SBATCH --ntasks=64

echo "let the games begin!"
for n in 1 2 4 8 16 32 64; do
	for x in 707 1000 1414 2000 2378 2828; do
		echo "== $n processes $x points =="
		mpirun -np "$n" "par_tri" "test"$x"pts93" | grep "Gent"
	done
done
