CXX = gcc

CXXFLAGS = -O3 -msse3 -ftree-vectorize -funroll-loops

FGDO_DIR = ..
EVALUATION_DIR = $(FGDO_DIR)/evaluation
SEARCH_DIR = $(FGDO_DIR)/searches
UTIL_DIR = $(FGDO_DIR)/util

APP_DIR = ../examples

SUM_OF_SQUARES = \
	$(APP_DIR)/sum_of_squares.o

EVALUATION_OBJS = \
	$(EVALUATION_DIR)/simple_evaluator.o \
	$(EVALUATION_DIR)/search_manager.o

SEARCH_OBJS = \
	$(SEARCH_DIR)/gradient.o \
	$(SEARCH_DIR)/hessian.o \
	$(SEARCH_DIR)/line_search.o \
	$(SEARCH_DIR)/newton_method.o \
	$(SEARCH_DIR)/outliers.o \
	$(SEARCH_DIR)/population.o \
	$(SEARCH_DIR)/recombination.o \
	$(SEARCH_DIR)/regression.o \
	$(SEARCH_DIR)/search_arguments.o \
	$(SEARCH_DIR)/search_log.o \
	$(SEARCH_DIR)/search_parameters.o \
	$(SEARCH_DIR)/synchronous_gradient_descent.o \
	$(SEARCH_DIR)/synchronous_newton_method.o
UTIL_OBJS = \
	$(UTIL_DIR)/io_util.o \
	$(UTIL_DIR)/matrix.o \
	$(UTIL_DIR)/settings.o

PROGS = sum_of_squares 

all: $(PROGS)

sum_of_squares: $(SUM_OF_SQUARES) $(SEARCH_OBJS) $(EVALUATION_OBJS) $(UTIL_OBJS)
	$(CXX) $(CXXFLAGS) -lm -Wall -o sum_of_squares $(SUM_OF_SQUARES) $(SEARCH_OBJS) $(EVALUATION_OBJS) $(UTIL_OBJS)

.c.o:
	$(CXX) $(CXXFLAGS) -Wall -x c -c $< -o $@

.C.o:
	$(CXX) $(CXXFLAGS) -Wall -x c -c $< -o $@
clean:
	rm -f $(PROGS) $(SUM_OF_SQUARES) $(EVALUATION_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) 
