/*****************************************************************************
Programmer: Alex Dyhdalo
Last Mod: 8/2017
*****************************************************************************/


#ifndef AUX_H
#define AUX_H

#include <iostream>
#include <string>
#include "matrix_define.h"
#include <cstdio>

// struct to hold the input file variables
struct parameters 
{
	size_t basis;
	size_t particles;
	double hw;
    std::string m_ref;
    std::string m_mat;
	bool N_flag;
	bool O_flag;
	bool P_flag;
  	bool Q_flag;
  	bool G_flag;
	bool diag_toggle;
	bool two_body_toggle;
	std::string sdp_file;
	std::string diag_file;
	bool mscheme_toggle;
};

// struct to hold the constraint flag toggles
struct con_flags
{
	bool N_flag;
	bool O_flag;
	bool P_flag;
  	bool Q_flag;
  	bool G_flag;
	bool diag_toggle;
	bool two_body_toggle;
};

parameters read_in_inputs ();
double kron_del(const size_t i, const size_t j);

template <typename Array> 
void print(std::ostream& os, const Array & A);

template <> void print<double> (std::ostream& os, const double & x);
template <> void print<int> (std::ostream& os, const int & x);

size_t Q_loop (size_t);
size_t G_loop (size_t);

void header_sdp_file (con_flags flag_pass, size_t, size_t Nnum, size_t Onum, size_t Pnum, size_t Qnum, size_t Gnum, std::ofstream & sdp_out);

double degen (size_t);

#endif
