unroll=$1
kernels=$(find autogen/* -name "*.c")
for kernel in ${kernels}; do
	kernel=$(basename ${kernel} ".c")
	echo "generating $kernel"
	./generate.sh $kernel ${unroll}
done
