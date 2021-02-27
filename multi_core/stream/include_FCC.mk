CC   = fcc
GCC  = fcc
LINKER = $(CC)

ifeq ($(ENABLE_OPENMP),true)
OPENMP   = -fopenmp
endif

CFLAGS   = -O3 -Ksimd_reg_size=512 -march=armv8.2-a+sve -Knounroll -Iincludes ${LIKWID_FLAG} $(OPENMP)
OPTIONS += -O3 -Ksimd_reg_size=512 -march=armv8.2-a+sve -Knounroll -Iincludes ${LIKWID_FLAG} $(OPENMP)
LFLAGS   = $(OPENMP)
DEFINES  = -D_GNU_SOURCE
INCLUDES =
LIBS     =
