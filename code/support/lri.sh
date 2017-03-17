#!/bin/sh
#SBATCH --mail-type = ALL
#SBATCH --mail-user = kogaluu@ccf.org
#SBATCH --job-name  = benchmark
#SBATCH --time      = 1-00:00:00

## Options are described below:
## -N 8 <request 8 nodes for the job>
## mpirun <execute the job using Open MPI> 
## -np 1 <execute only one copy of the program>
## --map_by ppr:1:node <launch processes one per node>
## --bind-to board <allow the slave to access all cores on the node> 
## --report-bindings <diagnostic trace of the binding>

salloc -N 8 mpirun -np 1 --map-by ppr:1:node --bind-to board --report-bindings R CMD BATCH --no-save combine.R
