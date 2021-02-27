CC   = gcc
GCC  = gcc
LINKER = $(CC)

ifeq ($(ENABLE_OPENMP),true)
OPENMP   = -fopenmp
endif

CFLAGS   = -Ofast -ffreestanding $(OPENMP)
CFLAGS  += -msve-vector-bits=512 -fno-unroll-loops -march=armv8.2-a+sve -funroll-loops
CFLAGS  +=  -Wl,-T/opt/FJSVxos/mmm/util/bss-2mb.lds -L/opt/FJSVxos/mmm/lib64 -lmpg 
LFLAGS   = $(OPENMP)
DEFINES  = -D_GNU_SOURCE
OPTIONS  += -msve-vector-bits=512 -fno-unroll-loops -march=armv8.2-a+sve -funroll-loops -Ofast -ffreestanding
INCLUDES =
LIBS     =
