mpirun -np 32 -machinefile machinefiles/mpi_machinefile ./mpi_astronomy -asynch -stars $1 -parameters $2 -nm -s $3 -nm_type line_search -nm_evaluations $4 -nm_iterations $5
