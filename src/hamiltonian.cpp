#include "hamiltonian.h"
#include <sstream>
#include <fstream>


/*
template<> void print<double>(std::ostream& os, const double & x)
{
  os << x << "\t";
}
*/


/***************************************************************


 Function to populate the 1-body part of the Hamiltonian. 
 The 1-body part here is T + U for kinetic energy T and external potential U


***************************************************************/


void populate_1body (const two_array & ref_m, two_array & h1_mat, const std::string obme_filename)
{

  size_t mat_length = h1_mat.size();

  int obme_size = get_obme_lines(obme_filename);

//  std::cout << "OBME SIZE " << obme_size << "\n";

  two_array obme (boost::extents[obme_size][6]);

  read_in_obme (obme, obme_filename);

  for (size_t i = 0; i < mat_length; ++i)
  {
    double n1  = ref_m [i][1];
    double l1  = ref_m [i][2];
    double j1  = ref_m [i][3];
    double mj1 = ref_m [i][4];

  for (size_t j = 0; j < mat_length; ++j)
  {
    double n2  = ref_m [j][1];
    double l2  = ref_m [j][2];
    double j2  = ref_m [j][3];
    double mj2 = ref_m [j][4];

    double matrix_element = 0.;

    if (l1 == l2 and j1 == j2 and mj1 == mj2) matrix_element = find_obme_me (n1, n2, l1, j1, mj1, obme);

	h1_mat[i][j] = matrix_element;
  }
  }

}

/***************************************************************



***************************************************************/

double find_obme_me (const double n1, const double n2, const double l, const double j, const double mj, two_array & obme)
{
  size_t mat_extent = obme.size();

  for (size_t i = 0; i < mat_extent; i++)
  {
    double obme_n1 = obme [i][0];
    double obme_n2 = obme [i][1];
    double obme_l  = obme [i][2];
    double obme_j  = obme [i][3];
    double obme_mj = obme [i][4];
    double obme_me = obme [i][5];


    if (n1 == obme_n1 and n2 == obme_n2 and l == obme_l and j == obme_j and mj == obme_mj)
      return obme_me;
  }


  return 0.;
}

/***************************************************************



***************************************************************/

int get_obme_lines (const std::string obme_filename)
{
  const char * reference_file = obme_filename.c_str();
  // input file stream for m_scheme
  std::ifstream ref_in (reference_file);
 
  size_t total_lines = 0;
  std::string dummy;

  // find total number of defined reference lines
  while (std::getline (ref_in, dummy)) ++total_lines;


  // clear the file stream, reset to read in the elements
  ref_in.clear();
  ref_in.seekg (0, std::ios::beg);

  return total_lines;
}


/***************************************************************


 Function to read in the data file for the one-body matrix elements.


***************************************************************/


void read_in_obme (two_array & obme, const std::string obme_filename)
{
  const char * m_reference_file = obme_filename.c_str();
  // input file stream for m_scheme
  std::ifstream ref_in (m_reference_file);
 
  size_t total_lines = 0;
  std::string dummy;
  double n1, n2, l, j, mj, me;

  // find total number of defined reference lines
  while (std::getline (ref_in, dummy))
  ++total_lines;


  // clear the file stream, reset to read in the elements
  ref_in.clear();
  ref_in.seekg (0, std::ios::beg);

  size_t ref_size = obme.size();
  size_t ele_in = 0;

  // read in and assign the references line by line
  // to the matrix ref_m
  for (size_t i = 0; i < total_lines; i++)
  {
  std::string orbit_dummy_1;
  std::string orbit_dummy_2;
  std::getline (ref_in, dummy);

  if (!dummy.length() || dummy[0] == '#')     // skip zero length lines and lines that start with #
    continue;

  std::stringstream ss;

  ss << dummy;            // read in the line to stringstream ss


  ss >> n1 >> n2 >> l >> j >> mj >> me;  // assign values of the line




    obme [ele_in][0] = n1;    // reference number of the line
    obme [ele_in][1] = n2;          // principle quantum number
    obme [ele_in][2] = l;          // orbital angular mom.
    obme [ele_in][3] = j;    // total angular mom.
    obme [ele_in][4] = mj;  // total angular mom. projection
    obme [ele_in][5] = me;         // isospin projection (should all be -1.0)
    ele_in++;

  if (ele_in >= ref_size)
    break;

  }
  

  return;
}


/***************************************************************

Wrapper function to read in single-particle m scheme reference 
file, populate 1-body Hamiltonian matrix elements, and then 
populate 2-body matrix elements.

***************************************************************/


void fullm_populate_hamiltonian (two_array & array_ref_obme, two_array & array_ref_tbme, two_array & h1_mat, two_array & h2_mat, const std::string ref_obme, const std::string me_obme, const std::string ref_tbme, const std::string me_tbme, const bool two_body_toggle, const int nmax) 
{

    read_in_reference_m_scheme (array_ref_obme, ref_obme);
    std::cout << "OBME REFERENCE READ" << std::endl;
    populate_1body (array_ref_obme, h1_mat, me_obme);
    std::cout << "1 BODY POPULATED" << std::endl;
    if (two_body_toggle)
    {
      std::cout << "TBME REFERENCE READ" << std::endl;
      populate_2body (array_ref_obme, array_ref_tbme, h2_mat, me_tbme, nmax);
      std::cout << "2 BODY POPULATED" << std::endl;
    }

}


/***************************************************************



***************************************************************/

size_t total_tbme_states (const std::string tbme_filename)
{
  const char * ref_file = tbme_filename.c_str();
  std::ifstream ref_in (ref_file);
 
  size_t total_lines = 0;
  std::string dummy;

  // find total number of defined reference lines
  while (std::getline (ref_in, dummy)) ++total_lines;

  return total_lines;
}

/***************************************************************



***************************************************************/

void readin_ref_tbme (two_array & ref_tbme, const std::string tbme_filename)
{

  const char * ref_file = tbme_filename.c_str();
  // input file stream for m_scheme
  std::ifstream ref_in (ref_file);
 
  size_t total_lines = 0;
  std::string dummy;
  double num, n1, l1, j1, mj1, n2, l2, j2, mj2, sp1, sp2;

  // find total number of defined reference lines
  while (std::getline (ref_in, dummy))
  ++total_lines;


  // clear the file stream, reset to read in the elements
  ref_in.clear();
  ref_in.seekg (0, std::ios::beg);

  size_t ref_size = ref_tbme.size();
  size_t ele_in = 0;

  // read in and assign the references line by line
  // to the matrix ref_m
  for (size_t i = 0; i < total_lines; i++)
  {
  	std::getline (ref_in, dummy);

  	if (!dummy.length() || dummy[0] == '#')     // skip zero length lines and lines that start with #
    	continue;

  	std::stringstream ss;

 	ss << dummy;            // read in the line to stringstream ss


  	ss >> num >> n1 >> l1 >> j1 >> mj1 >> n2 >> l2 >> j2 >> mj2 >> sp1 >> sp2;  // assign values of the line


    ref_tbme [ele_in][0] = num;    // reference number of the line

    ref_tbme [ele_in][1] = n1;          // principle quantum number
    ref_tbme [ele_in][2] = l1;          // orbital angular mom.
    ref_tbme [ele_in][3] = j1;    // total angular mom.
    ref_tbme [ele_in][4] = mj1;   // total angular mom. projection

    ref_tbme [ele_in][5] = n2;          // principle quantum number
    ref_tbme [ele_in][6] = l2;          // orbital angular mom.
    ref_tbme [ele_in][7] = j2;    // total angular mom.
    ref_tbme [ele_in][8] = mj2;   // total angular mom. projection
    ref_tbme [ele_in][9] = sp1;
    ref_tbme [ele_in][10] = sp2;

    ele_in++;

    if (ele_in >= ref_size)
    	break;

  }
  

  return;

}

/***************************************************************



***************************************************************/

void populate_2body (const two_array & ref_m, const two_array & ref_tbme, two_array & h2_mat, const std::string me_tbme, const int nmax)
{
  // input file stream
  const char * matrix_file = (me_tbme).c_str();
  std::ifstream matrix_in (matrix_file);

  size_t total_lines = 0;
  std::string dummy;
  int alpha, beta, gamma, delta;
  double value;

  // get the total number of lines in the matrix elements file
  while(std::getline(matrix_in, dummy))  ++total_lines;
 
  // clear and reset the file stream
  matrix_in.clear();
  matrix_in.seekg(0, std::ios::beg);

  const size_t obme_size = ref_m.size();
  const size_t tbme_size = h2_mat.size();

//  std::cout << "OBME SIZE " << obme_size << std::endl;
//  std::cout << "TBME SIZE " << tbme_size << std::endl;


  // set the size of the m_scheme holder
  //m_scheme_matrix.set_size (m_size, m_size);


  for (size_t n = 0; n < total_lines; n++)
  {
//  	  std::cout << "LOOP TOP" << std::endl;

	  std::getline (matrix_in, dummy);

	  if (!dummy.length() || dummy[0] == '#')
	    continue;


	  bool h1 = false;
	  bool h2 = false;
	  bool h3 = false;
	  bool h4 = false;

	  int n1 = -1;
	  int n2 = -1;
	  int n3 = -1;
	  int n4 = -1;

	  int l1 = -1;
	  int l2 = -1;
	  int l3 = -1;
	  int l4 = -1;

	  double j1 = -1;
	  double j2 = -1;
	  double j3 = -1;
	  double j4 = -1;

	  double mj1 = -1;
	  double mj2 = -1;
	  double mj3 = -1;
	  double mj4 = -1;


	  std::stringstream ss;

	  ss << dummy;

	  ss >> alpha >> beta >> gamma >> delta >> value;

	  // check to see if alpha, beta, gamma, delta we're reading in
	  // matches the reference numbers in our reference matrix
	  // if yes --> read it in, otherwise discard it

	  for (size_t w = 0; w < obme_size; w++)
	  {
	    if (alpha == ref_m[w][0]) 
	    {
	      h1  = true;
	      n1  = ref_m[w][1];
	      l1  = ref_m[w][2];
	      j1  = ref_m[w][3];
	      mj1 = ref_m[w][4];
	    }

	    if (beta == ref_m[w][0]) 
	    {
	      h2 = true;
	      n2  = ref_m[w][1];
	      l2  = ref_m[w][2];
	      j2  = ref_m[w][3];
	      mj2 = ref_m[w][4];
	    }

	    if (gamma == ref_m[w][0]) 
	    {
	      h3 = true;
	      n3  = ref_m[w][1];
	      l3  = ref_m[w][2];
	      j3  = ref_m[w][3];
	      mj3 = ref_m[w][4];
	    }

	    if (delta == ref_m[w][0]) 
	    {
	      h4 = true;
	      n4  = ref_m[w][1];
	      l4  = ref_m[w][2];
	      j4  = ref_m[w][3];
	      mj4 = ref_m[w][4];
	    }
	  }
	  
	  if (h1 and h2 and h3 and h4)
	  {
	  	int left = -10;
	  	int right = -10;

	  	if (nmax < (2*n1 + l1 + 2*n2 + l2) or nmax < (2*n3 + l3 + 2*n4 + l4))
	  	{
//	  	  std::cout << "HIT" << std::endl; 
	  	  continue;
	  	}

//	  	std::cout << "HIT" << std::endl;

	  	for (size_t loop = 0; loop < tbme_size; loop++)
	  	{
	  		int n1_loop  = ref_tbme[loop][1];
	  		int l1_loop  = ref_tbme[loop][2];
	  		double j1_loop  = ref_tbme[loop][3];
	  		double mj1_loop = ref_tbme[loop][4];

	  		int n2_loop  = ref_tbme[loop][5];
	  		int l2_loop  = ref_tbme[loop][6];
	  		double j2_loop  = ref_tbme[loop][7];
	  		double mj2_loop = ref_tbme[loop][8];

	  		if (n1 == n1_loop and n2 == n2_loop and l1 == l1_loop and l2 == l2_loop 
	  			and j1 == j1_loop and j2 == j2_loop and mj1 == mj1_loop and mj2 == mj2_loop)
	  			left = ref_tbme[loop][0] - 1;


	  		if (n3 == n1_loop and n4 == n2_loop and l3 == l1_loop and l4 == l2_loop 
	  			and j3 == j1_loop and j4 == j2_loop and mj3 == mj1_loop and mj4 == mj2_loop)
	  			right = ref_tbme[loop][0] - 1;
	  	}

//	  	if (left < 0 or right < 0) std::cerr << "ERROR: BAD TBME MATRIX INDEX " << left << " " << right << std::endl;

	  	if (left >= 0 and right >=0) h2_mat [left][right] = value;
	  }



  }

  return;

}

/***************************************************************



***************************************************************/


void compactify_h2 (const two_array & ref_m, two_array & comp_h2, five_array & h2_mat, std::ofstream & diag_out, const bool diag_toggle)
{
  const size_t bsize = h2_mat.size();

  double value;

  size_t alpha, beta, gamma, delta;

  if(diag_toggle)
  {
    diag_out << "Output of nonzero 2-body ME" << std::endl << std::endl;    
  }

  for (size_t i = 0; i < bsize; ++i)
  {
  for (size_t j = 0; j < bsize; ++j)
  {
    if (i >= j)
      continue;

  for (size_t k = 0; k < bsize; ++k)
  {
  for (size_t l = 0; l < bsize; ++l)
  {
    if (k >= l)
      continue;

            value = h2_mat [i][j][k][l][0];

            alpha = h2_mat [i][j][k][l][1];
            beta  = h2_mat [i][j][k][l][2];
            gamma = h2_mat [i][j][k][l][3];
            delta = h2_mat [i][j][k][l][4];

            if (diag_toggle)
            {
              if (value != 0.)
              {
              	diag_out << alpha << " " << beta << " " << gamma << " " << delta << " " << value << " \t";

                for (size_t loop = 0; loop < bsize; loop++)
                {
                  if (alpha == ref_m[loop][0])
                    diag_out << ref_m[loop][6] << "\t";
                }


                for (size_t loop = 0; loop < bsize; loop++)
                {
                  if (beta == ref_m[loop][0])
                    diag_out << ref_m[loop][6] << "\t";
                }


                for (size_t loop = 0; loop < bsize; loop++)
                {
                  if (gamma == ref_m[loop][0])
                    diag_out << ref_m[loop][6] << "\t";
                }


                for (size_t loop = 0; loop < bsize; loop++)
                {
                  if (delta == ref_m[loop][0])
                    diag_out << ref_m[loop][6] << "\n";
                }

              }
            }

            size_t ips = i + 1;
            size_t jps = j + 1;
            size_t kps = k + 1;
            size_t lps = l + 1;

            size_t left  = jps - ips + (2*bsize - ips) * (ips - 1)/2 - 1;
            size_t right = lps - kps + (2*bsize - kps) * (kps - 1)/2 - 1;


  /*         if (i == 0 && j == 4 && k == 0)
            {
            	std::cout << ref_m[l][1] << "\t" << ref_m[l][2] << "\n";
            }

           if (i == 0 && j == 1)
            {
            	std::cout << right << " " << l << " " << k << "\t" << value << "\n"; 
            }
*/

//            std::cout << value << std::endl;

            comp_h2 [left][right] = value;

  }
  }
  }
  }

  if (diag_toggle)
    diag_out << std::endl << std::endl;
}



/***************************************************************

 
 Function to read in the reference file for m scheme.


***************************************************************/


void read_in_reference_m_scheme (two_array & ref_m, const std::string m_ref_file)
{
  const char * m_reference_file = (m_ref_file).c_str();
  // input file stream for m_scheme
  std::ifstream m_ref_in (m_reference_file);
 
  if (m_ref_in.fail())
  throw "ERROR: Cannot open m-scheme reference file";
 
  size_t total_lines = 0;
  std::string dummy;
  double ref_num, n, l, j , m_j, tz;

  // find total number of defined reference lines
  while (std::getline (m_ref_in, dummy))
  ++total_lines;


  // clear the file stream, reset to read in the elements
  m_ref_in.clear();
  m_ref_in.seekg (0, std::ios::beg);

  size_t ref_size = ref_m.size();
  size_t ele_in = 0;

  // read in and assign the references line by line
  // to the matrix ref_m
  for (size_t i = 0; i < total_lines; i++)
  {
  std::string orbit_dummy_1;
  std::string orbit_dummy_2;
  std::getline (m_ref_in, dummy);

  if (!dummy.length() || dummy[0] == '#')     // skip zero length lines and lines that start with #
    continue;

  std::stringstream ss;

  ss << dummy;            // read in the line to stringstream ss

  // MORTEN VS. HEIKO READ IN FILE FORMAT

  ss >> orbit_dummy_1 >> orbit_dummy_2 >> ref_num >> n >> l >> j >> m_j >> tz;  // assign values of the line

//  ss >> ref_num >> n >> l >> j >> m_j >> tz;  // assign values of the line

//  std::cout << ref_num << " " << n << " " << l << " " << j << " " << m_j << " " << tz << "\n";


  // only extract neutron-neutron states
  if (tz == -1.)
  {
    ref_m [ele_in][0] = ref_num;    // reference number of the line
    ref_m [ele_in][1] = n;          // principle quantum number
    ref_m [ele_in][2] = l;          // orbital angular mom.
    ref_m [ele_in][3] = j * 0.5;    // total angular mom.
    ref_m [ele_in][4] = m_j * 0.5;  // total angular mom. projection
    ref_m [ele_in][5] = tz;         // isospin projection (should all be -1.0)
    ref_m [ele_in][6] = ele_in;     // new index for sp orbital
    ele_in++;
  }

  if (ele_in >= ref_size)
    break;

  }
  



  return;
}

/***************************************************************

Function to read in the actual m scheme matrix  elements 
calculated from Morten's code. 

THIS WILL BREAK IF THE INPUT FILE .dat CHANGES TITLE OR FORMAT.

***************************************************************/

void read_in_matrix_m_scheme (const two_array & ref_m, five_array & h2_mat, const std::string m_mat_file)
{
  // input file stream
  const char * m_matrix_file = (m_mat_file).c_str();
  std::ifstream m_matrix_in (m_matrix_file);

  size_t total_lines = 0;
  std::string dummy;
  int alpha, beta, gamma, delta;
  double value;

  if (m_matrix_in.fail())
  throw "ERROR: Cannot read m-scheme matrix file";

  // get the total number of lines in the matrix elements file
  while(std::getline(m_matrix_in, dummy))
  ++total_lines;
 
  // clear and reset the file stream
  m_matrix_in.clear();
  m_matrix_in.seekg(0, std::ios::beg);

  const size_t m_size = h2_mat.size();

  // set the size of the m_scheme holder
  //m_scheme_matrix.set_size (m_size, m_size);


  int i,j,k,l;

  for (size_t n = 0; n < total_lines; n++)
  {
  std::getline (m_matrix_in, dummy);

  if (!dummy.length() || dummy[0] == '#')
    continue;

  i = -1;
  j = -1;
  k = -1;
  l = -1;

  std::stringstream ss;

  ss << dummy;

  ss >> alpha >> beta >> gamma >> delta >> value;

  // check to see if alpha, beta, gamma, delta we're reading in
  // matches the reference numbers in our reference matrix
  // if yes --> read it in, otherwise discard it

  for (size_t w = 0; w < m_size; w++)
  {
    if (alpha == ref_m[w][0]) 
      i = w;

    if (beta == ref_m[w][0]) 
      j = w;

    if (gamma == ref_m[w][0]) 
      k = w;

    if (delta == ref_m[w][0]) 
      l = w;
  }
  

  if (i >= 0 && j >= 0 && k >= 0 && l >= 0) // only activates if all 4 "if" statements above are true
    {
      h2_mat [i][j][k][l][0] = value;
      h2_mat [i][j][k][l][1] = alpha;
      h2_mat [i][j][k][l][2] = beta;
      h2_mat [i][j][k][l][3] = gamma;
      h2_mat [i][j][k][l][4] = delta;
    }
  }

  // reset index term on reference matrix to span 0 to m_size instead of e.g. 3, 4, 9, 10, etc... for neutrons

//  for (size_t n = 0; n < m_size; n++)
//  ref_m[n][0] = n;

//  print(std::cout, ref_m);

  return;
}