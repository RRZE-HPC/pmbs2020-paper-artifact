export LC_NUMERIC="en_US.UTF-8"

nthreads=12
mem_blup=24
unroll=1
outFolder="multicore_scaling"
mkdir -p ${outFolder}
outFile="${outFolder}/u${unroll}.txt"
./generate.sh "stencil_2d5pt" ${unroll}
printf "%8s, %14s, %14s\n" "Threads" "MLUPs" "MB/s" > ${outFile}
for thread in `seq 1 ${nthreads}`; do
	#OMP_NUM_THREADS=$thread OMP_PLACES=cores OMP_PROC_BIND=close  ./bench_stencil_2d5pt ${unroll} 10000 > tmp.txt
    likwid-pin -q -c S0:0-$((thread-1)) ./bench_stencil_2d5pt ${unroll} 10000 > tmp.txt
	mlups=$(grep "MLUPS" tmp.txt | cut -d"=" -f2)
	mbs=$(echo "${mlups}*${mem_blup}" | bc -l)
	printf "%8d, %14.4f, %14.4f\n" ${thread} ${mlups} ${mbs} >> ${outFile}
done
