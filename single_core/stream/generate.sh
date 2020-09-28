kernel="$1"
inc_line=$(grep "${kernel}(" kernels/"$kernel.cpp")
num_streams=$(echo $inc_line | grep -o "," | wc -l)
kernel_call=$(grep -o "${kernel}(.*" kernels/"$kernel.cpp" | sed "s#double##g" | sed "s#*##g" | sed "s#int##g")

sed -e "s#@FWD_DECLARATION@#${inc_line};#g" main.cpp > generated_main.cpp
sed -i "s#@N_ARRAYS@#${num_streams}#g" generated_main.cpp
sed -i "s#@KERNEL@#${kernel_call};#g" generated_main.cpp

unroll_fac="$2"
cat kernels/${kernel}.cpp | sed -e "s@#UNROLL_FACTOR_SUB#@$unroll_fac@g"  > tmp_${kernel}.cpp
perl ghost_unroll.pl tmp_${kernel}.cpp > ${kernel}.cpp

./compile.sh "$kernel"

rm -rf tmp_${kernel}.cpp ${kernel}.cpp ${kernel}.o
rm -rf generated_main.cpp 
rm -rf ${kernel}.s
