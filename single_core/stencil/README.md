##How to compile and run
* For compiling with unrolling factor of 8  run: `./generate.sh stencil_2d5pt 8`
* Execute `./bench_stencil_2d5pt 8`
* To run specific sizes  use the optional third argument, i.e., `./bench_stencil_2d5pt 8 1000` will run
stencil of dimension 500x1000. Outer to inner dimension is in the ratio of 1:2.
* The data-set size, time and cycles are reported in the output.
