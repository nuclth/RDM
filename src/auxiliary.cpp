/*****************************************************************************
Programmer: Alex Dyhdalo
Last Mod: 8/2017
*****************************************************************************/

#include "auxiliary.h"
#include <sstream>
#include <fstream>

/**********************************************

Function which reads in the user specified values in
inputs.inp to an instance of the struct parameters.

**********************************************/

parameters read_in_inputs ()
{
  struct parameters input;				// local struct to hold user values
 
  std::ifstream input_in ("inputs.inp");		// the relevant input file stream
  std::string dummy;
  int counter = 0;					// counter lets us know which variable we're reading in

  while (std::getline (input_in, dummy)) 		// while runs as long as there are lines to read in the input file
  {
	if (!dummy.length() || dummy[0] == '#')		// skip zero length lines and lines that start with #
	  continue;

	std::stringstream ss;

	ss << dummy;					// read in the line to the stringstream ss

	if (counter == 0)				// read in parameters by the way they're ordered in "inputs.inp"
	  ss >> input.basis;

	if (counter == 1)
	  ss >> input.particles;

	if (counter == 2)
	  ss >> input.hw;

	if (counter == 3)
	  ss >> input.m_ref;

  if (counter == 4)
    ss >> input.m_mat;

  if (counter == 5)
    ss >> input.two_body_toggle;

  if (counter == 6)
    ss >> input.N_flag;

  if (counter == 7)
    ss >> input.O_flag;

  if (counter == 8)
    ss >> input.P_flag;

  if (counter == 9)
    ss >> input.Q_flag;

  if (counter == 10)
    ss >> input.G_flag;

  if (counter == 11)
    ss >> input.diag_toggle;  

  if (counter == 12)
    ss >> input.sdp_file;  

  if (counter == 13)
    ss >> input.diag_file;  

  if (counter == 14)
    ss >> input.mscheme_toggle;  

	counter++;					// increase the counter after reading in a parameter	
  }

  return input;						// return the struct 
}


/***************************************************************

Kronecker delta function. Returns 1 if i=j and 0 otherwise.

***************************************************************/

double kron_del(const size_t i, const size_t j)
{

  if (i == j)
    return 1.;

  return 0.;
}

/***************************************************************


 Function to print out an arbitrary multiarray


***************************************************************/

template <typename Array>
void print(std::ostream& os, const Array & A)
{
  typename Array::const_iterator i;
  os.precision(2);
  os << "[";
  for (i = A.begin(); i != A.end(); ++i) {
    print(os, *i);
//    if (boost::next(i) != A.end())
//      os << ',' << "\t";
  }
  os << "]";
  os << "\n";
}

/***************************************************************


 Individual print functions depending on whether template takes 
 double or int.


***************************************************************/


template<> void print<double>(std::ostream& os, const double & x)
{
  os << x << "\t";
}

template<> void print<int>(std::ostream& os, const int & x)
{
  os << x;
}

/***************************************************************

Quick and dirty loop function to figure out total number of Q 
constraints from basis size.

***************************************************************/

size_t Q_loop (size_t bsize)
{
  size_t Q_num = 0;

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = 0; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = 0; l < bsize; l++)    // loop over matrix column
  {
    Q_num++;
  }
  }
  }
  }

/*  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = i+1; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = i; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = k+1; l < bsize; l++)    // loop over matrix column
  {
    if (j < l && k <= i)
      continue;

    Q_num++;
  }
  }
  }
  }*/

  return Q_num;
}

/***************************************************************

Quick and dirty loop function to figure out total number of G 
constraints from basis size.

***************************************************************/

size_t G_loop (size_t bsize)
{
  size_t G_num = 0;

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = 0; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = 0; l < bsize; l++)    // loop over matrix column
  {

    G_num++;
  }
  }
  }
  }

/*  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = 0; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = j; l < bsize; l++)    // loop over matrix column
  {
    
    if (j == l && k < i)
      continue;

    G_num++;
  }
  }
  }
  }*/

  return G_num;
}


/***************************************************************

Function to output number of constraints, number of sub blocks,
and the sub block dimensions to the top of the SDP file.

***************************************************************/

void header_sdp_file (con_flags flag_pass, size_t bsize, size_t Nnum, size_t Onum, size_t Pnum, size_t Qnum, size_t Gnum, std::ofstream & sdp_out)
{

  size_t cons = 0;
  size_t blocks = 0;

  if (flag_pass.N_flag)
  {
    cons += Nnum;
    blocks++;
  }

  if (flag_pass.O_flag)
  {
    cons += Onum;
    blocks++;
  }

  if (flag_pass.two_body_toggle)
    blocks++;
  
  if (flag_pass.P_flag)
    cons += Pnum;

  if (flag_pass.Q_flag)
  {
    cons += Qnum;
    blocks++;
  }

  if (flag_pass.G_flag)
  {
    cons += Gnum;
    blocks++;
  }


  sdp_out << cons << std::endl;
  sdp_out << blocks << std::endl;

  if (flag_pass.N_flag)
  {
    sdp_out << bsize << " ";
  }

  if (flag_pass.O_flag)
  {
    sdp_out << bsize << " ";
  }
  
  if (flag_pass.two_body_toggle)
  {
    sdp_out << bsize * bsize << " ";
  }

  if (flag_pass.Q_flag)
  {
      sdp_out << bsize * bsize << " ";
//    sdp_out << bsize * (bsize-1)/2 << " ";
  }

  if (flag_pass.G_flag)
  {
    sdp_out << bsize * bsize << " ";
  }

  sdp_out << std::endl;

}

double degen (size_t num)
{
	double factor = 0;

	if (num == 1)
		factor = 2;

	else if (num == 2)
		factor = 4;

	else if (num == 3)
		factor = 2;

	else if (num == 4)
		factor = 6;

	else if (num == 5)
		factor = 4;		

	else if (num == 6)
		factor = 2;

	else if (num == 7)
		factor = 4;

	else if (num == 8)
		factor = 2;	

	else if (num == 9)
		factor = 6;	

	else if (num == 10)
		factor = 4;	

	return factor;
}
