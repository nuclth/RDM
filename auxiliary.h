#ifndef auxiliary_h
#define auxiliary_h

#include <cstdio>
#include <string>
#include "arraydefs.h"

struct parameters 
{
  size_t basis;
  size_t particles;
  double hw;
  std::string m_ref;
  std::string m_mat;
};

struct con_flags
{
  bool N_flag;
  bool O_flag;
  bool P_flag;
  bool Q_flag;
  bool G_flag;
  bool two_body_toggle;
  bool diag_toggle;
  bool redundant_check;
};

parameters read_in_inputs ();

void diag_format (std::ofstream & diag_out, const two_array & c_matrix, const size_t N_num, const three_array & N_con, const one_array & N_val, const size_t O_num, const three_array & O_con, const one_array & O_val, const size_t P_num, const three_array & P_con, const one_array & P_val, const size_t Q_num, const three_array & Q_con, const one_array & Q_val, const size_t G_num, const three_array & G_con, const one_array & G_val, const con_flags flag_pass);

template <typename Array>
void print(std::ostream& os, const Array & A);

template<> void print<double>(std::ostream& os, const double & x);

template<> void print<int>(std::ostream& os, const int & x);

size_t full_matrix_extent (const size_t bsize, const con_flags flag_pass);

#endif