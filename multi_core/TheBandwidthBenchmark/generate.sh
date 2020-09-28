kernel="$1"
unroll_fac="$2"
cat autogen/${kernel}.c | sed -e "s@#UNROLL_FACTOR_SUB#@$unroll_fac@g"  > tmp_${kernel}.c
perl ghost_unroll.pl tmp_${kernel}.c > src/${kernel}.c

rm -rf tmp_${kernel}.cpp 
