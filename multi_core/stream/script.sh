nthreads=48
make

out_folder="size_4GB"
mkdir -p ${out_folder}
out_file="${out_folder}/u8.txt"
for (( threads=1; threads<=${nthreads}; threads++ )); do
	OMP_PLACES=cores OMP_PROC_BIND=close OMP_NUM_THREADS=$threads taskset -c 0-$((threads-1)) ./bwbench-GCC > tmp.txt
	bw=$(cat tmp.txt | tail -n +7 | head -n -2 | cut -d":" -f2 | cut -d"," -f1 | tr "\n" ",")
	bench=$(cat tmp.txt | tail -n +7 | head -n -2 | cut -d":" -f1 | tr "\n" ",")
	if [[ "${threads}" == "1" ]]; then
		echo "#Threads,${bench}" > ${out_file}
	fi
	echo "${threads},${bw}" >> ${out_file}
done
