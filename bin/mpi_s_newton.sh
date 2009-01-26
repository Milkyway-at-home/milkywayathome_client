mpirun -np 32 -machinefile machinefiles/mpi_machinefile ./mpi_astronomy -synch -stars stars-$1.txt -parameters parameters-$1.txt -nm -nm_iterations $2
