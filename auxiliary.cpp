
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

	if (counter == 0)
	  ss >> input.basis;

	if (counter == 1)
	  ss >> input.particles;

	if (counter == 2)
	  ss >> input.hw;

  	if (counter == 3)
      ss >> input.m_ref;

  	if (counter == 4)
      ss >> input.m_mat;


	counter++;					// increase the counter after reading in a parameter	
  }

  return input;						// return the struct 
}

/***************************************************************

 Function to print out an arbitrary array

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

 Function to output the constraint matrices to output.

***************************************************************/


void diag_format (std::ofstream & diag_out, const two_array & c_matrix, const size_t N_num, const three_array & N_con, const one_array & N_val, const size_t O_num, const three_array & O_con, const one_array & O_val, const size_t P_num, const three_array & P_con, const one_array & P_val, const size_t Q_num, const three_array & Q_con, const one_array & Q_val, const size_t G_num, const three_array & G_con, const one_array & G_val, const con_flags flag_pass)
{



    diag_out << "H Constraint Matrix" << std::endl << std::endl;

    print(diag_out, c_matrix);

    diag_out << std::endl << std::endl;


    if (flag_pass.N_flag)
    {
      diag_out << N_num << " terms - " << "N constraint matrix output" << std::endl << std::endl;

      print(diag_out, N_con);
      
      diag_out << std::endl << std::endl;



      diag_out << "N constraint values" << std::endl << std::endl;

      print(diag_out, N_val);

      diag_out << std::endl << std::endl;
    }


    if (flag_pass.O_flag)
    {
      diag_out << O_num << " terms - " << "O constraint matrix output" << std::endl << std::endl;

      print(diag_out, O_con);
      
      diag_out << std::endl << std::endl;



      diag_out << "O constraint values" << std::endl << std::endl;

      print(diag_out, O_val);

      diag_out << std::endl << std::endl;
    }


    if (flag_pass.P_flag)
    {
      diag_out << P_num << " terms - " << "P constraint matrix output" << std::endl << std::endl;

      print(diag_out, P_con);
      
      diag_out << std::endl << std::endl;



      diag_out << "P constraint values" << std::endl << std::endl;

      print(diag_out, P_val);

      diag_out << std::endl << std::endl;
    }


  

    if (flag_pass.Q_flag)
    {
      diag_out << Q_num << " terms - " << "Q constraint matrix output" << std::endl << std::endl;

      print(diag_out, Q_con);
      
      diag_out << std::endl << std::endl;



      diag_out << "Q constraint values" << std::endl << std::endl;

      print(diag_out, Q_val);

      diag_out << std::endl << std::endl;
    }




    if (flag_pass.G_flag)
    {
      diag_out << G_num << " terms - " << "G constraint matrix output" << std::endl << std::endl;

      print(diag_out, G_con);
      
      diag_out << std::endl << std::endl;



      diag_out << "G constraint values" << std::endl << std::endl;

      print(diag_out, G_val);

      diag_out << std::endl << std::endl;
    }
 
}


/***************************************************************

 Function to determine the size of each matrix in the problem.

***************************************************************/

size_t full_matrix_extent (const size_t bsize, const con_flags flag_pass)
{

  size_t cmat = 2 * bsize;

  if (flag_pass.two_body_toggle)
  {
      cmat = 2*bsize + bsize * (bsize-1)/2;

    if (flag_pass.Q_flag)
      cmat = 2*bsize + 2 * bsize * (bsize-1)/2;

    if (flag_pass.G_flag)
      cmat = 2*bsize + 2 * bsize * (bsize-1)/2 + bsize * bsize;
  }

  return cmat;
}