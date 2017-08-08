#include "jscheme.h"


/***************************************************************

 
 Function to read in the reference file for m scheme.


***************************************************************/


void read_in_reference_j_scheme (two_array & ref_j, const std::string j_ref_file, std::ofstream & diag_out, const bool diag_toggle)
{
  const char * j_reference_file = (j_ref_file).c_str();
  // input file stream for m_scheme
  std::ifstream j_ref_in (j_reference_file);
 
  if (j_ref_in.fail())
  throw "ERROR: Cannot open m-scheme reference file";
 
  size_t total_lines = 0;
  std::string dummy;
  double ref_num, n, l, j , m_j, tz;

  // find total number of defined reference lines
  while (std::getline (j_ref_in, dummy))
  ++total_lines;


  // clear the file stream, reset to read in the elements
  j_ref_in.clear();
  j_ref_in.seekg (0, std::ios::beg);

  size_t ref_size = ref_j.size();
  size_t ele_in = 0;

  // read in and assign the references line by line
  // to the matrix ref_m
  for (size_t i = 0; i < total_lines; i++)
  {
  std::string orbit_dummy_1;
  std::string orbit_dummy_2;
  std::getline (j_ref_in, dummy);

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
    ref_j [ele_in][0] = ref_num;    // reference number of the line
    ref_j [ele_in][1] = n;          // principle quantum number
    ref_j [ele_in][2] = l;          // orbital angular mom.
    ref_j [ele_in][3] = j * 0.5;    // total angular mom.
    ref_j [ele_in][4] = m_j * 0.5;  // total angular mom. projection
    ref_j [ele_in][5] = tz;         // isospin projection (should all be -1.0)
    ref_j [ele_in][6] = ele_in;     // new index for sp orbital
    ele_in++;
  }

  if (ele_in >= ref_size)
    break;

  }
  
  
  
  if(diag_toggle)
  {
    diag_out << "Single particle orbitals pulled from REF file" << std::endl << std::endl;
//    print(diag_out, ref_m);            // print the resulting matrix
    diag_out << std::endl << std::endl;
  }


  return;
}

/***************************************************************

Function to read in the actual m scheme matrix  elements 
calculated from Morten's code. 

THIS WILL BREAK IF THE INPUT FILE .dat CHANGES TITLE OR FORMAT.

***************************************************************/

void read_in_matrix_j_scheme (const two_array & ref_j, five_array & h2_mat, const std::string j_mat_file)
{
  // input file stream
  const char * j_matrix_file = (j_mat_file).c_str();
  std::ifstream j_matrix_in (j_matrix_file);

  size_t total_lines = 0;
  std::string dummy;
  int alpha, beta, gamma, delta;
  double value;

  if (j_matrix_in.fail())
  throw "ERROR: Cannot read m-scheme matrix file";

  // get the total number of lines in the matrix elements file
  while(std::getline(j_matrix_in, dummy))
  ++total_lines;
 
  // clear and reset the file stream
  j_matrix_in.clear();
  j_matrix_in.seekg(0, std::ios::beg);

  const size_t j_size = h2_mat.size();

  // set the size of the m_scheme holder
  //m_scheme_matrix.set_size (m_size, m_size);


  int i,j,k,l;

  for (size_t n = 0; n < total_lines; n++)
  {
  std::getline (j_matrix_in, dummy);

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

  for (size_t w = 0; w < j_size; w++)
  {
    if (alpha == ref_j[w][0]) 
      i = w;

    if (beta == ref_j[w][0]) 
      j = w;

    if (gamma == ref_j[w][0]) 
      k = w;

    if (delta == ref_j[w][0]) 
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