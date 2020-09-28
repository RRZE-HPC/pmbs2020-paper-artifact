# Artifact Description: Performance Modeling of Streaming Kernels and Sparse Matrix-Vector Multiplication on A64FX

* Christie Alappat, RRZE University of Erlangen-Nürnberg, christie.alappat@fau.de
* Thomas Gruber, RRZE University of Erlangen-Nürnberg, thomas.gruber@fau.de
* Georg Hager, RRZE University of Erlangen-Nürnberg, georg.hager@fau.de
* Jan Laukemann, RRZE University of Erlangen-Nürnberg, jan.laukemann@fau.de
* Nils Meyer, University of Regensburg, nils.meyer@ur.de
* Gerhard Wellein, RRZE University of Erlangen-Nürnberg, gerhard.wellein@fau.de
* Tilo Wettig, University of Regensburg, tilo.wettig@ur.de

## A.1 Abstract
The A64FX CPU powers the current #1 supercomputer in the Top500 list.
Although it is a traditional cache-based multicore processor, its peak performance and 
memory bandwidth rival accelerator devices.
Code for such a new architecture requires a good grasp of its performance features as
a basis for performance models.
In this paper we construct the Execution-Cache-Memory (ECM) performance model for the
A64FX processor in the FX700 supercomputer and validate it using streaming benchmarks.
Applying the model to sparse matrix-vector multiplication (SpMV), we motivate why the
CSR matrix storage format is inappropriate and how the SELL-C-σ format with suitable
code optimizations can achieve bandwidth saturation with SpMV.


## A.2 Description
### A.2.1 Check-list (artifact meta information)
- Compilation: gcc
- Binary: aarch64
- Hardware: A64FX (FX700)
- Publicly available?: yes

### A.2.2 How software can be obtained (if available)
Check out this repository.
```
git clone https://github.com/RRZE-HPC/pmbs2020-paper-artifact
```

### A.2.3 Hardware dependencies
We ran on a FX700 system (A64FX architecture) at 1.8 GHz with 48 cores per CPU.
The results should be reproducible on any FX700 system.

### A.2.4 Software dependencies
* GCC 10.1.1

### A.2.5 Datasets
None necessary, everything is part of the code.


## A.3 Installation
None necessary.


## A.4 Experiment workflow
To validate our results use the following commands.

Download script and benchmark codes:
```
git clone https://github.com/RRZE-HPC/pmbs2020-paper-artifact
cd pmbs20-paper-appendix/
```

**TODO...**


## A.5 Evaluation and expected result
**TODO in-core**
Compare numbers to Table II.

**TODO ECM**
Compare numbers to Table III

**TODO full node**
Compare numbers to Figure 5 (right).


## A.6 Experiment customization
None


## A.7 Notes
None
