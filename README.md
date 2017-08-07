# Reduced Density Matrices (RDMs) for Neutron Drops using Semi-Definite Programming (SDP)

A program to output an SDP data format file for neutrons in a harmonic trap interacting via minnesota potential. This can then be solved using an external SDP solver to find the ground state energy (along with the associated one- and two-body RDMs). Only m scheme choice of basis is currently supported.

## Details 

Source files are kept in the ```src``` directory, input m scheme matrix elements in the ```me_files``` directory, and output data files in the ```sdp_files``` directory.

The SDP data file is outputted in sparse data format as indicated by the 

```-s```

at the end of the file. Note that the sparse format (and solvers that we use in general) assumes that all constraint matrices are symmetric. Output format is described below. The term in the () indicates the relevant number in the included **sdp_example.dat-s** file.

* The total number of constraints
* The number of sub-blocks 
* The relevant sizes of each sub-block in order
* The values of each of the constraint matrices
* The output of each matrix with only the non-zero upper triangle specified
	* which constraint matrix is being specified
	* which block is specified
	* the relevant row
	* the relevant column
	* the value of the matrix 


## Instructions

* Compile and link the program by running the command ```make``` 

* Change any inputs as desired in **inputs.inp**

* Run with command ```./run_rdm```

* SDP runs currently available for two solvers: csdp or sdpa. For an output data file with name **output_file.dat-s** run this using

```csdp output_file.dat-s```

or

```sdpa output_file.dat-s solution.dat```

where solution.dat will store the solution to the SDP. For other options, try ```man sdpa``` or refer to relevant manual files.

## Dependencies



## Issues


