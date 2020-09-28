LIKWID_INC ?= -I/usr/local/include
LIKWID_DEFINES ?= -DLIKWID_PERFMON
LIKWID_LIB ?= -L/usr/local/lib

ifeq ($(strip $(ENABLE_LIKWID)),true)
INCLUDES += -I${LIKWID_INCDIR}
DEFINES +=  -DLIKWID_PERFMON
LIBS += -llikwid
LFLAGS += -L${LIKWID_LIBDIR}
endif
