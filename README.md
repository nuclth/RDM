# Reduced Density Matrices (RDMs) for Neutron Drops using Semi-Definite Programming (SDP)

A program to output an SDP data format file for neutrons in a harmonic trap interacting via minnesota potential. This can then be solved using an external SDP solver to find the ground state energy (along with the associated one- and two-body RDMs). Only m scheme choice of basis is currently supported.

## Details 

Source files are kept in the ```src``` directory, input m scheme matrix elements in the ```me_files``` directory, and output data files in the ```sdp_files``` directory.

The SDP data file is outputted in sparse data format as indicated by the 

```-s```

at the end of the file. Note that the sparse format (and solvers that we use in general) assumes that all constraint matrices are symmetric and sub-block formats between matrices are identical. Output format is given as:

* The total number of constraints
* The number of sub-blocks 
* The relevant sizes of each sub-block in order (pure consectuive diagonal entires are given with a - prefix)
* The values of each of the constraint matrices
* A line for each non-zero entry in the upper triangle part of a given constraint matrix
	* which constraint matrix is being specified starting from zero
	* which block is specified
	* the relevant row
	* the relevant column
	* the value of the matrix 

Note that rows and columns reset for each block and start at 1. For example, supposing the matrix

```
[[6 1 0 0 0]
 [1 8 0 0 0]
 [0 0 4 0 0]
 [0 0 0 3 0]
 [0 0 0 0 7]]
```

was our first and only constraint matrix term with value 4.5, the SDP file would be given by,

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

## Instructions

* Compile and link the program by running the command ```make``` 

* Change any inputs as desired in **inputs.inp**

* Run with command ```./run_rdm```

* SDP runs currently available for two solvers: csdp or sdpa. For an output data file with name **output_file.dat-s** run this using

```csdp output_file.dat-s```

or

```sdpa output_file.dat-s solution.dat```

where solution.dat will store the solution to the SDP. For other options, try ```man sdpa``` or refer to relevant manual files.

## Branches

### Dynamic

* master - main development branch

* jscheme - development branch to allow for use of a J scheme basis. Should allow for much larger basis sets when working. Open issues about anti-symmetry, normalization of matrix elements, construction of constraint matrices, etc...

### static

* blockdiag - Branch used for project to put SDP file into local block diagonal structure for m scheme. The idea was that the potential for different ```Jz``` and/or different parities (say ```Jz = 1``` and ```Jz = 0``` both with even parity) do not talk to one another. That is, all of these matrix elements are zero. This should speed up the SDP solver immensely. However, these goals were not achieved. The energy behaved erratically and I was unable to reproduce the results of the antisym part of the code. My working theory on why this is the following: the ground state of the system is in good total ```J``` but by rotational invariance is independent of total ```Jz``` (may be wrong). Therefore, different ```Jz``` states can contribute to the ground state and by making them block diagonal, we are cutting off this coupling in the system (even if individual matrix elements do not couple). I am not convinced by this explanation (or at least by the way it is stated now). However, I still believe the idea here is sound and should be applied to a J scheme calculation when up and running. The system has good total ```J``` and parity which can be extracted easily, so the blockdiagonal structure should work there. 
	* codeclean - Offshoot of blockdiag. First attempt to clean up my code (use different headers/source files, more commenting, etc...). Lots of futile work here as this branch was ultimately abadoned once blockdiag continued to not work.

* antisym - First branch to take advantage of the anti-symmetric nature of the 2RDM and potential matrix elements. Allows for smaller matrix sizes with the same basis. Quite a good amount of speed up against the original implementation. For a basis size of r, the original implementation used matrices of rank ```r*r``` while the antisym exploit allows for rank ```r*(r-1)/2``` 

* original_imp - Points to the point on master branch where the first 'working' implemenation was made in m scheme. Naive implementation with no benefits due to block diagonal structure, anti-symmetry of the potential, code maintenance, etc... Marker kept for purposes of comparison/bench marking with future work.

## Dependencies

N/A 

## Issues

Doesn't work yet.

## Benchmarks

Calculations done here with parameters,

```
hb*c = 197.326 MeV-fm
m    = 938.92 MeV
hb*w = 10 MeV
b    = hb*c / sqr(hb*w*m) = 2.0364 fm
```

using the **test_index.dat** and **test_mscheme_20.dat** data files (from Heiko's code).

Exact diagonlization of the m scheme Hamiltonian gives the energy:


```
2  states - n = 0, l = 0 - E = 25.5266
8  states - n = 0, l = 1 - E = 25.2015
18 states - n = 0, l = 2 - E = 25.1561
20 states - n = 1, l = 0 - E = 24.9853
26 states - n = 1, l = 1 - E = 24.927
36 states - n = 1, l = 2 - E = 24.9188 
38 states - n = 2, l = 0 - E = 24.9043
44 states - n = 2, l = 1 - E = 24.8995
54 states - n = 2, l = 2 - E = 24.8984
```

while the SDP formulation with P condition gives:

```
2  states - n = 0, l = 0 - E = 25.52655
8  states - n = 0, l = 1 - E = 25.20149 
18 states - n = 0, l = 2 - E = 25.15606 
20 states - n = 1, l = 0 - E = 24.82686
26 states - n = 0, l = 3 - E = 24.71367
36 states - n = 1, l = 1 - E = 24.69817
38 states - n = 1, l = 2 - E = 24.66184
44 states - n = 2, l = 0 - E = 
54 states - n = 2, l = 0 - E = 
```

The SDP formulation with PQ condition gives:

```
2  states - n = 0, l = 0 - E = 
8  states - n = 0, l = 1 - E = 
18 states - n = 0, l = 2 - E = 
20 states - n = 1, l = 0 - E = 24.82686
26 states - n = 0, l = 3 - E = 
36 states - n = 1, l = 1 - E = 
38 states - n = 1, l = 2 - E = 
44 states - n = 2, l = 0 - E = 
54 states - n = 2, l = 0 - E = 
```

## Benchmarks OLD

For 2 particles interacting via Minnesota potential (with usual parameters, **add paper ref**) and choice of:

```
hb*c = 197.326 MeV-fm
m    = 938.92 MeV
hb*w = 10 MeV
b    = hb*c / sqr(hb*w*m) = 2.0364 fm
```

the ground state energy of the system is ```E = 24.265 MeV```. Building up a basis with only closed shells using Morten's matrix elements, the energy progression is:

```
2  states - n = 0, l = 0 - E = 25.526540
8  states - n = 0, l = 1 - E = 25.201479
18 states - n = 0, l = 2 - E = 25.156046
20 states - n = 1, l = 0 - E = 24.826845
32 states - n = 0, l = 3 - E = 24.826768 REDO
40 states - n = 1, l = 1 - E = 24.727879
50 states - n = 1, l = 2 - E = 24.727862
52 states - n = 2, l = 0 - E = 24.727862
```

while Heiko's matrix elements give,

```
2  states - n = 0, l = 0 - E = 25.526560
8  states - n = 0, l = 1 - E = 25.201498
18 states - n = 0, l = 2 - E = 25.156066 
20 states - n = 1, l = 0 - E = 24.826269
34 states - n = 0, l = 3 - E = 24.823097
40 states - n = 1, l = 1 - E = 24.710750
50 states - n = 1, l = 2 - E = 24.695374
52 states - n = 2, l = 0 - E = 24.659119
```

**UPDATE**: Morten's matrix elements are given above for the triangular (Nmax truncation) model space. Going to a larger space or using square (Emax truncation) model space reproduces Heiko's numbers.

More benchmarks using 2 matrix full trace condition (NN). 2 particles in 8 basis states (Heiko code):

```
P        - E = 25.201496
P NN     - E = 25.201496
P NN Q   - E = 25.201496
P NN Q G - E = 25.201496

NN       - E = 23.027203
NN Q     - E = 25.198266
NN Q G   - E = 25.201496
```
