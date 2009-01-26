mpirun -np 32 -machinefile machinefiles/mpi_machinefile ./mpi_astronomy -synch -stars $1 -parameters $2 -gd -gd_min_threshold 0.00001 -gd_iterations $3
