#ifndef AUXILIARY_H
#define AUXILIARY_H

#include <iostream>
#include <string>

// struct to hold the input file variables
struct parameters 
{
	size_t basis;
	size_t particles;
    std::string ref_obme;
    std::string me_obme;
    std::string ref_tbme;
    std::string me_tbme;
};

// struct to hold the constraint flag toggles
struct con_flags
{
	bool N_flag;
	bool O_flag;
	bool P_flag;
  	bool Q_flag;
  	bool G_flag;

  	bool T1_flag;
  	bool T2_flag;

  	bool NN_flag;

	bool diag_toggle;
	bool two_body_toggle;
};

parameters read_in_inputs ();

inline double kron_del(const size_t i, const size_t j)
{
  if (i == j)
    return 1.;

  return 0.;
}

template <typename Array> void print(std::ostream& os, const Array & A);

template <> void print (std::ostream& os, const int & x);
size_t T2_count (const size_t bsize);
size_t T2_DIM_count (const size_t bsize);

#endif
