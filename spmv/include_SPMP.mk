SPMP_INCDIR ?= ./SpMP-master
SPMP_DEFINES ?= -DUSE_SPMP
SPMP_LIB ?= ./SpMP-master/libspmp.a 

ifeq ($(strip $(ENABLE_SPMP)),true)
INCLUDES += -I${SPMP_INCDIR}
DEFINES += ${SPMP_DEFINES}
LIBS += ${SPMP_LIB}
endif
