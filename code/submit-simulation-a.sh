#!/bin/bash

for SEPARATION in 1 2 3 4 5 6 7 8 9 10
do
	sbatch submit-simulation-a ${SEPARATION}
done

echo "Simulation setting A. All jobs submitted."
