CC=gcc
FLAGS=-std=c99
LIB=-lm -lmkl_sequential -lmkl_core -lmkl_intel_lp64 -liomp5
EXECUTABLE=meta_vc_est
SOURCES=param_est.c

all: meta_vc_est

meta_vc_est: $(SOURCES)
	$(CC) $(FLAGS) $(LIB) $(SOURCES) -o $(EXECUTABLE)
