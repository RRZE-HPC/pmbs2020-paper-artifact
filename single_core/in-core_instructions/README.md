## How to compile and run
* Run `make` on a A64FX CPU in this directory. The .so files should be build in `./SVE/`.
* Execute `./ibench SVE/ 1.8`
* The results for each instruction form is printed to stdout.  
  Each instruction form name is encoded in the pattern "`mnemonic-operands-[LAT|TP]`".  
  * `TP` stands for "throughput", while `LAT` stands for "latency" of the instruction form.
  * `d`, `p`, `zd`, `x` stand for the according register types.
  * `m` stands for memory address and is followed by `b` for "base register" and optionally `i` for "base register + index"
