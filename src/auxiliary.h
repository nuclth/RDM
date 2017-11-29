#ifndef AUXILIARY_H
#define AUXILIARY_H

#include <iostream>
#include <string>

using std::string;

// struct to hold the input file variables
struct parameters 
{
	size_t basis_hw;
	size_t nmax;
	size_t particles;

	bool N_flag;
	bool O_flag;
	bool P_flag;

	bool two_body_toggle;
};

// struct 
struct string_holder
{

  string morten_spm;

  string ref_obme;
  string me_obme;

  string ref_tbme;
  string me_tbme;

  string no_flag;
  string h2_flag;
  string pflag_info;
};

// struct to hold the constraint flag toggles
struct con_flags
{
	bool N_flag;
	bool O_flag;
	bool P_flag;

	bool two_body_toggle;
};

parameters read_in_inputs ();

string_holder string_reader (const size_t nmax, const size_t basis_hw);

inline double kron_del(const size_t i, const size_t j)
{
  if (i == j)
    return 1.;

  return 0.;
}



#endif
