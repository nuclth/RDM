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
	  ss >> input.ref_obme;

 	if (counter == 3)
      ss >> input.me_obme;

  	if (counter == 4)
      ss >> input.ref_tbme;

  	if (counter == 5)
  	  ss >> input.me_tbme;

	counter++;					// increase the counter after reading in a parameter	
  }

  return input;						// return the struct 
}


size_t T2_count (const size_t bsize)
{
	size_t T2_num = 0;

  	for (size_t i = 0;   i < bsize; i++)      
 	{
 	for (size_t j = 0;   j < bsize; j++)     
 	{
  	for (size_t k = j+1; k < bsize; k++)   
  	{
  	for (size_t l = 0;   l < bsize; l++)   
  	{
  	for (size_t m = 0;   m < bsize; m++)   
  	{
  	for (size_t n = m+1; n < bsize; n++)    
  	{
  		T2_num++;
  	}
  	}
	}
	}
	}
	}

	return T2_num;
}


size_t T2_DIM_count (const size_t bsize)
{
	size_t dim = 0;

	    for (size_t ip = 0;    ip < bsize; ip++)      // loop over ith constraint matrix
	    {
	    for (size_t jp = 0;    jp < bsize; jp++)      // loop over jth constraint matrix
	    {
	    for (size_t kp = jp+1; kp < bsize; kp++)    // loop over matrix row
	    {

	    	dim++;
	    }
		}
		}

	return dim;
}
/***************************************************************

Kronecker delta function. Returns 1 if i=j and 0 otherwise.

***************************************************************/

/*inline double kron_del(const size_t i, const size_t j)
{

  if (i == j)
    return 1.;

  return 0.;
}*/

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



template<> void print<int>(std::ostream& os, const int & x)
{
  os << x;
}

