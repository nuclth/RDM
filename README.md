# Reduced Density Matrices (RDMs) for Neutron Drops using Semi-Definite Programming (SDP)

A program to output an SDP data format file for neutrons in a harmonic trap interacting via the Minnesota potential. This can then be solved using an external SDP solver to find the ground state energy (along with the associated one- and two-body RDMs). Only the m scheme harmonic oscillator basis is currently supported.

## Usage Instructions

There are three steps that must be followed before running the RDM program (which creates the SDP file).

1. Create one-body matrix elements using the Mathematica notebook in the ```src``` directory. Simply open the notebook and run it. Note however that the output filename, basis $\hbar \omega$, $n_{\text{max}}$, and $l_{\text{max}}$ currently all needs to be input by hand. Choices of $n_{\text{max}}$ and $l_{\text{max}}$ should match some choice of total $N$ for $N_{\text{max}}$ truncation, $N = 2n + l$. 

2. Create two-body matrix elements using Morten's ME code and then put the TBME output file and single-particle output file in the correct directories (```me_files/tbme/``` and ```me_files/ref_files/``` respectively). Again, the filenames, basis $\hbar \omega$, $n_{\text{max}}$, and $l_{\text{max}}$ all need to be put in by hand. These are all specified in the ```renorm.ini``` file. Then run the program using ```./vrenorm.exe```. 

3. Create the single-particle orbital states list, two-body basis states list, and data files with SDP constraint information using the python file ```nmax_count.py```. The user needs to specify the ```nmax``` variable in the file (where the ```nmax``` variable is our $N_{\text{max}}$ truncation) and then run the python file.

After this is all done, the user specifies the $N_{\text{max}}$ truncation, basis $\hbar \omega$, and the number of particles in the ```inputs.inp``` file (constraint flags in the inputs file are not yet operational; changing these currently needs to be done in the ```top_rdm.cpp``` file). Then make the program with ```make``` and run with ```./run_rdm```.

After the SDP file is generated (currently only with filename ```sdp_files/test_sdp.dat-s```), the SDP solver can then be called on it using for example:

* ```csdp test_sdp.dat-s```
* ```sdpa test_sdp.dat-s solution_out.dat```

## File Hierarchy and Description

Program source files are kept in the ```src``` directory, input m scheme matrix elements and single-particle files and two-particle basis files in the ```me_files``` directory, input flag files from the python program in the ```flag_files``` directory, and output data files in the ```sdp_files``` directory.

## SDP Data File Format

The SDP data file is outputted in sparse data format as indicated by the 

```-s```

at the end of the file. Note that the sparse format (and solvers that we use in general) assumes that all constraint matrices are symmetric and sub-block formats between matrices are identical. Output format is given as:

* The total number of constraints
* The number of sub-blocks 
* The relevant sizes of each sub-block in order (pure consecutive diagonal entires are given with a - prefix)
* The constant term for each constraint matrix relation
* A line for each non-zero entry in the upper triangle part of a given constraint matrix
	* which constraint matrix is being specified starting from zero for the matrix in the objective function
	* which matrix block is specified (starting from 1)
	* the relevant row (starting from 1)
	* the relevant column (starting from 1)
	* the value of the matrix 

Note that rows and columns reset for each block. For example, supposing the matrix

```
[[6 1 0 0 0]
 [1 8 0 0 0]
 [0 0 4 0 0]
 [0 0 0 3 0]
 [0 0 0 0 7]]
```

was our first and only constraint matrix term with the constant value 4.5 in the matrix relation. Our SDP file would be given by,

```
1
2
2 -3
4.5
0 1 1 1 6
0 1 1 2 1
0 1 2 2 8
0 2 1 1 4
0 2 2 2 3
0 2 3 3 7
```


## Benchmarks 

Calculations for benchmarks assume the following parameter values:


* $\hbar c = 197.326 \text{MeV-fm}$
* $m    = 938.92 \ \text{MeV}$
* $\hbar \omega = 10 \ \text{MeV}$
* $b    = \frac{\hbar c}{\sqrt{\hbar \omega m}} = 2.0364 \ \text{fm}$

