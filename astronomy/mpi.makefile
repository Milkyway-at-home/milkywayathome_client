BOINC_DIR = /projects/wcl/software/boinc

CXX = mpicc

MPI_ASTRONOMY = mpi_astronomy.o
SEPARATION = separation.o
ERRORS = errors.o
ASTRONOMY_OBJS = atSurveyGeometry.o numericalIntegration.o parameters.o probability.o stCoords.o stCnum.o stMath.o stVector.o star_points.o evaluation.o
SEARCH_OBJS = ../searches/gradient_descent.o ../searches/line_search.o ../searches/genetic_search.o ../searches/newton_method.o ../searches/differential_evolution.o ../searches/particle_swarm.o ../searches/synchronous_search.o ../searches/recombination.o ../searches/gradient.o ../searches/hessian.o ../searches/population.o
UTIL_OBJS = ../util/io_util.o ../util/matrix.o
EVALUATION_OBJS = ../evaluation/mpi_evaluator.o

PROGS = mpi_astronomy

all: $(PROGS)

mpi_astronomy: $(MPI_ASTRONOMY) $(ASTRONOMY_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) $(EVALUATION_OBJS)
	$(CXX) -lm -Wall -o mpi_astronomy $(MPI_ASTRONOMY) $(ASTRONOMY_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) $(EVALUATION_OBJS)

separation: $(SEPARATION) $(ASTRONOMY_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) $(EVALUATION_OBJS)
	$(CXX) -lm -Wall -o separation $(SEPARATION) $(ASTRONOMY_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) $(EVALUATION_OBJS)


errors: $(ERRORS) $(ASTRONOMY_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) $(EVALUATION_OBJS) 
	$(CXX) -lm -Wall -o errors $(ERRORS) $(ASTRONOMY_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) $(EVALUATION_OBJS)

.c.o:
	$(CXX) $(CXXFLAGS) -Wall -x c -c $< -o $@

.C.o:
	$(CXX) $(CXXFLAGS) -Wall -x c -c $< -o $@

clean:
	rm $(ASTRONOMY_OBJS) $(MPI_ASTRONOMY) $(SEARCH_OBJS) $(UTIL_OBJS) $(EVALUATION_OBJS)
