mpirun -np 32 -machinefile machinefiles/mpi_machinefile ./mpi_astronomy -asynch -stars stars-$1.txt -parameters parameters-$1.txt -nm -s $2 -nm_type line_search -nm_evaluations $3 -nm_iterations $4
