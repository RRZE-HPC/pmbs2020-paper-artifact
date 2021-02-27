# Supported: GCC, CLANG, ICC
TAG ?= GCC
ENABLE_OPENMP ?= true
ENABLE_LIKWID ?= true

#Feature options
OPTIONS  =  -DSIZE=130000000ull
OPTIONS +=  -DNTIMES=20
OPTIONS +=  -DARRAY_ALIGNMENT=64
OPTIONS +=  -ffast-math
OPTIONS +=  -DACLE_VERSION
#-fvariable-expansion-in-unroller --param max-variable-expansions-in-unroller=4
#OPTIONS +=  -DVERBOSE_AFFINITY
#OPTIONS +=  -DVERBOSE_DATASIZE
#OPTIONS +=  -DVERBOSE_TIMER
