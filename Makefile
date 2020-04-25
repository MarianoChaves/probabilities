# Makefile

# Name of the project
PROJ_NAME=main
 
# .cpp files
C_SOURCE=$(wildcard *.cpp)
 
# .h files
H_SOURCE=$(wildcard *.h)
 
# Object files
OBJ=$(C_SOURCE:.cpp=.o)
 
# Compiler
CC=g++
 
# Flags for compiler
prefix = /usr/local
exec_prefix = /usr/local
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

local_CFLAGS = -g -O4 -std=c++11

INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs)
preopen_modules:=$(shell $(globesconf) --dlpreopen)
local_LTLDFLAGS:=$(shell $(globesconf) --ltlibs)
local_GSLFLAGS:=$(gsl-config --libs)
local_GSFLAGS:=$(gsl-config --cflags)
ifdef preopen_modules
predefs = -DGLB_STATIC
endif

CC_FLAGS= -c -W -Wall -ansi -pedantic -g -O4 -std=c++11 
OMP_FLAGS= -fopenmp
GSL_FLAGS= -lgsl -lgslcblas
 

# Compilation and linking

all: $(PROJ_NAME)
 
$(PROJ_NAME): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS) $(local_GSFLAGS) $(local_LDFLAGS) $(local_GSLFLAGS) $(OMP_FLAGS) $(GSL_FLAGS)
 
%.o: %.cpp %.h
 
	$(CC) -o $@ $< $(CC_FLAGS) $(LDFLAGS) $(local_GSFLAGS) $(local_LDFLAGS) $(local_GSLFLAGS) $(INCFLAGS) $(OMP_FLAGS) $(GSL_FLAGS)
 
main.o: main.cpp $(H_SOURCE)
	$(CC) -o $@ $< $(CC_FLAGS) $(LDFLAGS) $(local_GSFLAGS) $(local_LDFLAGS) $(local_GSLFLAGS) $(INCFLAGS) $(OMP_FLAGS) $(GSL_FLAGS)
 
#clean:
#	 rm -rf *.o $(PROJ_NAME) *~
