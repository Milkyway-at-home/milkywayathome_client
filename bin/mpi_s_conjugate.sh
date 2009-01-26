mpirun -np 32 -machinefile machinefiles/mpi_machinefile ./mpi_astronomy -synch -stars $1 -parameters $2 -cgd -gd_reset 8 -gd_min_threshold 0.00001 -gd_iterations $3
