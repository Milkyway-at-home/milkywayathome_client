#ifndef GEM_MPI_EVALUATOR_H
#define GEM_MPI_EVALUATOR_H

void mpi_evaluator__init(int *argc, char*** argv);
void mpi_evaluator__read_data(void (*read_data)(int, int));
void mpi_evaluator__start();

int get_mpi_rank();

#endif
