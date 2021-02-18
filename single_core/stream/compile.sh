kernel="$1"
lkwid="off"
compiler="g++"
#compiler="FCC"
LIKWID_FLAG=""
LIKWID_LIB=""
if [[ "$likwid" == "on" ]]; then
	LIKWID_FLAG="-I${LIKWID_INCDIR} -DLIKWID_PERFMON"
	LIKWID_LIB="-L${LIKWID_LIBDIR} -llikwid"
fi
# g++
INCLUDE_FLAGS="-DTHPALLOC -DALIGNED -O3 -msve-vector-bits=512 -march=armv8.2-a+sve -Iincludes ${LIKWID_FLAG}"
# FCC
#INCLUDE_FLAGS="-DALIGNED -O3 -Ksimd_reg_size=512 -march=armv8.2-a+sve -Iincludes ${LIKWID_FLAG}"

${compiler} ${INCLUDE_FLAGS} -S ${kernel}.cpp
if [[ ${kernel} == "load" || ${kernel} == "load2" || ${kernel} == "load4" ]]; then
	cat ${kernel}.s > tmp_${kernel}.s
	cat tmp_${kernel}.s | sed -e "s@faddv.*@@g" | sed -e "s@fadd.*@@g" > ${kernel}.s
fi

${compiler} ${INCLUDE_FLAGS} -c ${kernel}.s
${compiler} ${INCLUDE_FLAGS} -c generated_main.cpp 
${compiler} ${INCLUDE_FLAGS} -c allocate.cpp 
${compiler} ${INCLUDE_FLAGS} ${kernel}.o allocate.o generated_main.o ${LIKWID_FLAG} ${LIKWID_LIB} -o bench_${kernel}
