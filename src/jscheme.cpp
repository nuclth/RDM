/*****************************************************************************
Programmer: Alex Dyhdalo
Last Mod: 8/2017
*****************************************************************************/

#include "jscheme.h"
#include "auxiliary.h"


/***************************************************************

 Function to read in the reference file for m scheme.

***************************************************************/


void read_in_reference_j_scheme (two_array & ref_j, const std::string j_ref_file, std::ofstream & diag_out, const bool diag_toggle)
{
  const char * j_reference_file = (j_ref_file).c_str();
  // input file stream for j_scheme
  std::ifstream j_ref_in (j_reference_file);
 
  if (j_ref_in.fail())
  throw "ERROR: Cannot open j-scheme reference file";
 
  size_t total_lines = 0;
  std::string dummy;
  double ref_num, n, l, twoj, tz, twonpl, HO_energy;

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
  std::getline (j_ref_in, dummy);

  if (!dummy.length() || dummy[0] == '#')     // skip zero length lines and lines that start with #
    continue;

  std::stringstream ss;

  ss << dummy;            // read in the line to stringstream ss

  // MORTEN VS. HEIKO READ IN FILE FORMAT

  ss >> orbit_dummy_1 >> ref_num >> n >> l >> twoj >> tz >> twonpl >> HO_energy;  // assign values of the line

//  ss >> ref_num >> n >> l >> j >> m_j >> tz;  // assign values of the line

//  std::cout << ref_num << " " << n << " " << l << " " << j << " " << m_j << " " << tz << "\n";


  // only extract neutron-neutron states
  if (tz == -1.)
  {
    ref_j [ele_in][0] = ref_num;    // reference number of the line
    ref_j [ele_in][1] = n;          // principal quantum number
    ref_j [ele_in][2] = l;          // orbital angular mom.
    ref_j [ele_in][3] = twoj * 0.5; // total angular mom.
    ref_j [ele_in][4] = tz;         // isospin projection (should all be -1.0)
    ref_j [ele_in][5] = twonpl;    // value of 2*n+l
    ref_j [ele_in][6] = HO_energy;  // HO energy for the sp state
    ele_in++;
  }

  if (ele_in >= ref_size)
    break;

  }
  
  
  
  if(diag_toggle)
  {
    diag_out << "Single particle orbitals pulled from REF file" << std::endl << std::endl;
//    print(diag_out, ref_j);            // print the resulting matrix
    diag_out << std::endl << std::endl;
  }


  return;
}

/***************************************************************

Function to read in the actual m scheme matrix elements 
calculated from Morten's code. 

THIS WILL BREAK IF THE INPUT FILE .dat CHANGES TITLE OR FORMAT.

***************************************************************/

void read_in_matrix_j_scheme (two_array & h2_mat, const two_array & twop_basis, const std::string j_mat_file)
{
  // input file stream
  const char * j_matrix_file = (j_mat_file).c_str();
  std::ifstream j_matrix_in (j_matrix_file);

  size_t total_lines = 0;
  std::string dummy;
  int tz, parity, twoJ;
  int alpha, beta, gamma, delta;
  double value;

  if (j_matrix_in.fail())
  throw "ERROR: Cannot read j-scheme matrix file";

  // get the total number of lines in the matrix elements file
  while(std::getline(j_matrix_in, dummy))
  ++total_lines;
 
  // clear and reset the file stream
  j_matrix_in.clear();
  j_matrix_in.seekg(0, std::ios::beg);

  const size_t j_size = h2_mat.size();

  // set the size of the m_scheme holder
  //m_scheme_matrix.set_size (m_size, m_size);


  int right, left;

  for (size_t n = 0; n < total_lines; n++)
  {
  std::getline (j_matrix_in, dummy);

  if (!dummy.length() || dummy[0] == '#')
    continue;

  right = -1;
  left = -1;

  std::stringstream ss;

  ss << dummy;

  ss >> tz >> parity >> twoJ >> alpha >> beta >> gamma >> delta >> value;

  // check to see if alpha, beta, gamma, delta we're reading in
  // matches the reference numbers in our twop basis matrix
  // if yes --> read it in, otherwise discard it

  for (size_t w = 0; w < j_size; w++)
  {
    if (alpha == twop_basis[w][0] and beta == twop_basis[w][1] and twoJ == twop_basis[w][2])
    {
      left = w;
    }

    if (alpha == twop_basis[w][1] and beta == twop_basis[w][0] and twoJ == twop_basis[w][2])
    {
      left = w;
    }

    if (gamma == twop_basis[w][0] and delta == twop_basis[w][1] and twoJ == twop_basis[w][2])
    {
      right = w;
    }

    if (gamma == twop_basis[w][1] and delta == twop_basis[w][0] and twoJ == twop_basis[w][2])
    {
      right = w;
    }
  }
  

  if (left >= 0 and right >= 0)  h2_mat [left][right] = value; // assign hamiltonian value

  }

  // reset index term on reference matrix to span 0 to m_size instead of e.g. 3, 4, 9, 10, etc... for neutrons

//  for (size_t n = 0; n < m_size; n++)
//  ref_m[n][0] = n;

//  print(std::cout, ref_m);

  return;
}

/***************************************************************

Function to read in J scheme single particle states and count the
total number of two-particle states for a given single particle
basis size.

***************************************************************/

size_t count_twopart_jscheme (const two_array & ref_j)
{
  size_t jcount = 0;
  size_t bsize = ref_j.size();

  for (size_t loop1 = 0;     loop1 < bsize; loop1++)
  {
  for (size_t loop2 = loop1; loop2 < bsize; loop2++)
  {

    int tz1 = ref_j [loop1][4];
    int tz2 = ref_j [loop2][4];

    if (tz1 != -1 or tz2 != -1) continue;

    size_t twoj1 = 2 * ref_j [loop1][3];
    size_t twoj2 = 2 * ref_j [loop2][3];

    size_t twojmax = twoj1 + twoj2;
    size_t twojmin = abs (twoj1 - twoj2);

//    std::cout << twojmin << " " << twojmax << std::endl;

    for (size_t loop3 = twojmin; loop3 <= twojmax; loop3+=2)
      jcount++; // end loop3

  } // end loop2
  } // end loop1

  return jcount;
}

/***************************************************************

Function to create an array list for all two particle basis
states in J scheme. Array has three entries:

1. ref number of single particle state 1
2. ref number of single particle state 2
3. total J value for sp state combo

***************************************************************/

void create_2pbasis_jscheme (two_array & twop_basis, const two_array & ref_j)
{
  size_t count = 0;

  size_t bsize = ref_j.size();

  for (size_t loop1 = 0;     loop1 < bsize; loop1++)
  {
  for (size_t loop2 = loop1; loop2 < bsize; loop2++)
  {

    int tz1 = ref_j [loop1][4];
    int tz2 = ref_j [loop2][4];

    if (tz1 != -1 or tz2 != -1) continue;

    size_t twoj1 = 2 * ref_j [loop1][3];
    size_t twoj2 = 2 * ref_j [loop2][3];

    size_t twojmax = twoj1 + twoj2;
    size_t twojmin = abs (twoj1 - twoj2);

//    std::cout << twojmin << " " << twojmax << std::endl;

    for (size_t twojtot = twojmin; twojtot <= twojmax; twojtot+=2)
    {
      twop_basis [count][0] = ref_j [loop1][0];
      twop_basis [count][1] = ref_j [loop2][0];
      twop_basis [count][2] = twojtot;
      count++;

    } // end twojtot loop

  } // end loop2
  } // end loop1

  return;
}

