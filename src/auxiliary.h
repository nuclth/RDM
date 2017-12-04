#ifndef AUXILIARY_H
#define AUXILIARY_H

#include <iostream>
#include <string>

using std::string;

struct parameters     // struct to hold the input file variables
{
	size_t basis_hw;
	size_t nmax;
	size_t particles;

	bool N_flag;
	bool O_flag;
	bool P_flag;

	bool two_body_toggle;
};

struct string_holder  // struct to hold strings for filenames
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
 

// function prototypes

parameters read_in_inputs ();
string_holder string_reader (const size_t nmax, const size_t basis_hw);
size_t total_obme_states (const string obme_filename);
size_t total_tbme_states (const string tbme_filename);
size_t count_no_flags (const string no_flag);
size_t count_P_cons (const string p_filename);

void init_con_blocks (const parameters flag_pass, const size_t NO_blocks, const size_t O_num, const size_t P_num, FILE * sdpa_out);


inline double kron_del(const size_t i, const size_t j)
{
  if (i == j)
    return 1.;

  return 0.;
}



#endif
