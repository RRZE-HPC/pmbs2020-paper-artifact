CC   = gcc
CXX  = g++
LINKER = $(CXX)

ifeq ($(ENABLE_OPENMP),true)
OPENMP   = -fopenmp
endif

CPPFLAGS   = -std=c++11 -Ofast -mcpu=native -fopenmp-simd -funroll-loops -fno-builtin -msve-vector-bits=512 -march=armv8.2-a+sve ${OPENMP} -Wno-write-strings
#-Ofast -ffreestanding -std=c99 $(OPENMP)
LFLAGS   = $(OPENMP)
DEFINES  = -D_GNU_SOURCE
INCLUDES =
LIBS     =

ENABLE_SCC ?= false
