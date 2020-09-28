## How to compile
* Configure using config.mk. 
* For using RCM: Intel SpMP library has to be installed, and `ENABLE_SPMP` has to be set to true in config.mk, along with paths in `include_SPMP.mk`
* For using LIKWID: LIKWID has to be installed, `ENABLE_LIKWID` has to bet to true in config.mk. Proper paths has to be set in `include_LIKWID.mk`
* After configuration run `make` to compiler


## How to run
* Get options using `./spmv-[TAG] -h`
* For running with HPCG matrix of size 128^3 on 1 CMG using SELL-32-64:
```
OMP_NUM_THREADS=$thread OMP_PLACES=cores OMP_PROC_BIND=close ./spmv-[TAG] -S 128 -C 32 -Z 64
```
* For other matrices, matrices can be read from MatrixMarket file using `-m` option.