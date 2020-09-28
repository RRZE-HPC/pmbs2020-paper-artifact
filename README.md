# Artifact Description: Performance Modeling of Streaming Kernels and Sparse Matrix-Vector Multiplication on A64FX

* Christie Alappat, RRZE University of Erlangen-Nürnberg, christie.alappat@fau.de
* Jan Laukemann, RRZE University of Erlangen-Nürnberg, jan.laukemann@fau.de
* Thomas Gruber, RRZE University of Erlangen-Nürnberg, thomas.gruber@fau.de
* Georg Hager, RRZE University of Erlangen-Nürnberg, georg.hager@fau.de
* Gerhard Wellein, RRZE University of Erlangen-Nürnberg, gerhard.wellein@fau.de
* Nils Meyer, University of Regensburg, nils.meyer@ur.de
* Tilo Wettig, University of Regensburg, tilo.wettig@ur.de

## A.1 Abstract
The A64FX CPU powers the current #1 supercomputer on the Top500
list.  Although it is a traditional cache-based multicore processor,
its peak performance and memory bandwidth rival accelerator devices.
Generating efficient code for such a new architecture requires a good understanding of its
performance features. Using these features, we construct the
Execution-Cache-Memory (ECM) performance
model for the A64FX processor in the FX700 supercomputer and
validate it using streaming loops. We also identify architectural peculiarities
and derive optimization hints. Applying the ECM model to sparse
matrix-vector multiplication (SpMV), we motivate why the CRS matrix
storage format is inappropriate and how the SELL-C-σ format with
suitable code optimizations can achieve bandwidth
saturation for SpMV. 


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
For SpMV we use matrices from SuiteSparse Matrix Collection
(https://suitesparse-collection-website.herokuapp.com).

## A.3 Installation
None necessary.


## A.4 Experiment workflow
To validate our results use the following commands.

Download script and benchmark codes:
```
git clone https://github.com/RRZE-HPC/pmbs2020-paper-artifact
cd pmbs20-paper-appendix/
```
The description for compiling and running can be found in the corresponding README.md files of each benchmark.


## A.5 Evaluation and expected result
### In-core instructions
Navigate to the [in-core_instructions directory](single_core/in-core_instructions):
```
cd single_core/in-core_instructions/
```
Run the validation accordingly to the README file and compare the results to Table II in the paper.

### Single-core
We provide code to reproduce the streaming kernels and the 2d five-point benchmark.  
Navigate to either one of the directories:
```
cd single_core/[stream|stencil]/
```
Run the validation accordingly to the README file and compare the results to Table III in the paper.

### Multi-core
We provide code to reproduce the streaming kernels and the 2d five-point benchmark.  
Navigate to either one of the directories:
```
cd multi_core/[stream|stencil]/
```
Run the validation accordingly to the README file and compare the results to Figure 4 in the paper.


## A.6 Experiment customization
None


## A.7 Notes
None
