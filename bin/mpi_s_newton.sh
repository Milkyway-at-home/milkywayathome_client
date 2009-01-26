mpirun -np 32 -machinefile machinefiles/mpi_machinefile ./mpi_astronomy -synch -stars $1 -parameters $2 -nm -nm_iterations $3
