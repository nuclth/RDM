#include "hamiltonian.h"
#include <sstream>
#include <fstream>
#include <map>


/***************************************************************

Wrapper function to read in single-particle m scheme reference 
file, populate 1-body Hamiltonian matrix elements, and then 
populate 2-body matrix elements.

***************************************************************/


void populate_hamiltonian (two_array & array_ref_obme, two_array & array_ref_tbme, two_array & h1_mat, two_array & h2_mat, const std::string ref_obme, const std::string me_obme, const std::string ref_tbme, const std::string me_tbme, const bool two_body_toggle, const int nmax, const std::string h2_flag) 
{

    readin_ref_obme (array_ref_obme, ref_obme);
    std::cout << "OBME REFERENCE READ" << std::endl;

    populate_1body (array_ref_obme, h1_mat, me_obme, nmax);
    std::cout << "1 BODY POPULATED" << std::endl;

    if (two_body_toggle)
    {
    	readin_ref_tbme (array_ref_tbme, ref_tbme);
    	std::cout << "TBME REFERENCE READ" << std::endl;

    	populate_2body (array_ref_obme, array_ref_tbme, h2_mat, me_tbme, nmax, h2_flag);
    	std::cout << "2 BODY POPULATED" << std::endl;
    }

}



/***************************************************************

 Function to read in the single-particle orbital reference file
 and populate array_ref_obme.

 ref_obme = "me_files/ref_files/nmax*_python_sp.dat"

***************************************************************/


void readin_ref_obme (two_array & array_ref_obme, const std::string ref_obme)
{
  const char * ref_file = (ref_obme).c_str();
  std::ifstream ref_in (ref_file);
 
  size_t total_lines = 0;
  std::string dummy;
  double num, n, l, j, mj, block;

  // find total number of defined reference lines
  while (std::getline (ref_in, dummy)) ++total_lines;


  // clear the file stream, reset to read in the elements
  ref_in.clear();
  ref_in.seekg (0, std::ios::beg);

  size_t ref_size = array_ref_obme.size();
  size_t ele_in = 0;

  // read in line by line and assign to array
  for (size_t i = 0; i < total_lines; i++)
  {

  	std::getline (ref_in, dummy);

 	// skip zero length lines and lines that start with #
	if (!dummy.length() || dummy[0] == '#')    
	    continue;

	std::stringstream ss;

 	ss << dummy;            
 	ss >> num >> n >> l >> j >> mj >> block;

    array_ref_obme [ele_in][0] = num;      // sp orbital number
    array_ref_obme [ele_in][1] = n;        // principal quantum number
    array_ref_obme [ele_in][2] = l;        // orbital angular mom.
    array_ref_obme [ele_in][3] = j * 0.5;  // total angular mom.
    array_ref_obme [ele_in][4] = mj * 0.5; // total angular mom. projection
    array_ref_obme [ele_in][5] = block;    // block number (non-zero for last spo in block, 0 otherwise)
    array_ref_obme [ele_in][6] = ele_in;   // index of array
    ele_in++;


	if (ele_in >= ref_size)
	    break;

  }
  

  return;
}


/***************************************************************


 Function to populate the 1-body part of the Hamiltonian, h1_mat. 
 The 1-body part here is T + U for kinetic energy T and external potential U.

 Arguments:

	array_ref_obme - matrix that holds info about the sp orbitals

	h1_mat - matrix to hold one-body matrix elments

 	obme_filename - data file with matrix elements of the one-body hamiltonian
 		e.g., "me_files/obme/nmax*_obme_hw*.dat"

	nmax - Nmax truncation value

***************************************************************/


void populate_1body (const two_array & array_ref_obme, two_array & h1_mat, const std::string obme_filename, const int nmax)
{

  const char * ref_file = obme_filename.c_str();
  std::ifstream ref_in (ref_file);
 
  size_t total_lines = 0;
  std::string dummy;
  int n1, n2, l;
  double j, mj, me;

  // find total number of defined reference lines
  while (std::getline (ref_in, dummy)) ++total_lines;

  // clear the file stream, reset to read in the elements
  ref_in.clear();
  ref_in.seekg (0, std::ios::beg);

  size_t obme_size = h1_mat.size();


  for (size_t i = 0; i < total_lines; i++)
  {
	std::getline (ref_in, dummy);

  	if (!dummy.length() || dummy[0] == '#')     
    	continue;

  	std::stringstream ss;

  	ss << dummy;            
  	ss >> n1 >> n2 >> l >> j >> mj >> me;  

  	if (nmax < n1 + l or nmax < n2 + l)
  		continue;

  	for (size_t loop1 = 0; loop1 < obme_size; loop1++)
  	{
  		int n1_ref     = array_ref_obme [loop1][1];
  		int l1_ref     = array_ref_obme [loop1][2];
  		double j1_ref  = array_ref_obme [loop1][3];
  		double mj1_ref = array_ref_obme [loop1][4];
  		int ele1       = array_ref_obme [loop1][6];

  	for (size_t loop2 = 0; loop2 < obme_size; loop2++)
  	{
  		int n2_ref     = array_ref_obme [loop2][1];
  		int l2_ref     = array_ref_obme [loop2][2];
  		double j2_ref  = array_ref_obme [loop2][3];
  		double mj2_ref = array_ref_obme [loop2][4];
  		int ele2       = array_ref_obme [loop2][6];

  		if (l1_ref != l2_ref or j1_ref != j2_ref or mj1_ref != mj2_ref)
  			continue;

  		if (n1 == n1_ref and n2 == n2_ref and l == l1_ref and j == j1_ref and mj == mj1_ref) h1_mat [ele1][ele2] = me;

  		if (n2 == n1_ref and n1 == n2_ref and l == l1_ref and j == j1_ref and mj == mj1_ref) h1_mat [ele2][ele1] = me;

  	} // end loop2
  	} // end loop1

  } // end outer for loop

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
  double num, n1, l1, n2, l2, sp1, sp2, block;

  double j1, mj1, j2, mj2;

  // find total number of defined reference lines
  while (std::getline (ref_in, dummy)) ++total_lines;


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


  	ss >> num >> n1 >> l1 >> j1 >> mj1 >> n2 >> l2 >> j2 >> mj2 >> sp1 >> sp2 >> block;  // assign values of the line


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
    ref_tbme [ele_in][11] = block;
    ref_tbme [ele_in][12] = ele_in;

    ele_in++;

    if (ele_in >= ref_size)
    	break;

  }
  

  return;

}

/***************************************************************



***************************************************************/

void populate_2body (const two_array & ref_m, const two_array & ref_tbme, two_array & h2_mat, const std::string me_tbme, const int nmax, const std::string h2_flag)
{

  const char * flag_file = (h2_flag).c_str();
  std::ifstream flag_in (flag_file);

  std::string dummy;
  size_t total_lines = 0;

  while (std::getline (flag_in, dummy)) ++total_lines;

  flag_in.clear();
  flag_in.seekg(0, std::ios::beg);

  two_array relational_sp (boost::extents[total_lines][10]);

  size_t b1, b2, m1, m2, sp1, sp2, sp3, sp4, msp1, msp2, msp3, msp4, block;

  for (size_t loop = 0; loop < total_lines; loop++)
  {

    std::getline (flag_in, dummy);

    if (!dummy.length() || dummy[0] == '#')
      continue;

    std::stringstream ss;

    ss << dummy;

    ss >> b1 >> b2 >> m1 >> m2 >> sp1 >> sp2 >> sp3 >> sp4 >> msp1 >> msp2 >> msp3 >> msp4 >> block;

    relational_sp [loop][0] = m1;
    relational_sp [loop][1] = m2;

    relational_sp [loop][2] = sp1;
    relational_sp [loop][3] = sp2;
    relational_sp [loop][4] = sp3;
    relational_sp [loop][5] = sp4;

    relational_sp [loop][6] = msp1;
    relational_sp [loop][7] = msp2;
    relational_sp [loop][8] = msp3;
    relational_sp [loop][9] = msp4;
  }

  total_lines = 0;


  const char * matrix_file = (me_tbme).c_str();
  std::ifstream matrix_in (matrix_file);



  // get the total number of lines in the matrix elements file
  while(std::getline(matrix_in, dummy))  ++total_lines;
 
  // clear and reset the file stream
  matrix_in.clear();
  matrix_in.seekg(0, std::ios::beg);

  size_t relational_size = relational_sp.size();

  std::map<int, int> morten_map;
  std::map<int, int> sign_map;

  for (size_t loop = 0; loop < relational_size; loop++)
  {
  	size_t msp1 = relational_sp [loop][6];
  	size_t msp2 = relational_sp [loop][7];

  	int m1 = relational_sp [loop][0];

  	int pair_map1 = cantor (msp1, msp2);
  	int pair_map2 = cantor (msp2, msp1);

  	morten_map.insert(std::make_pair(pair_map1, m1));
  	morten_map.insert(std::make_pair(pair_map2, m1));

  	sign_map.insert(std::make_pair(pair_map1,  1));
  	sign_map.insert(std::make_pair(pair_map2, -1));
  }

//  int cant_test1 = cantor (1, 2);
//  int cant_test2 = cantor (2, 6);
//  int cant_test3 = cantor (1, 38);

//  std::cout << morten_map.find(cant_test1)->second << std::endl;
//  std::cout << morten_map.find(cant_test2) << std::endl;
//  std::cout << morten_map.find(cant_test3) << std::endl;


//  two_array hold (boost::extents [total_lines][5]);

  for (size_t loop2 = 0; loop2 < total_lines; loop2++)
  {
//      std::cout << "LOOP TOP" << std::endl;

    size_t alpha, beta, gamma, delta;
    double value;

    std::getline (matrix_in, dummy);

    if (!dummy.length() || dummy[0] == '#')
      continue;



    std::stringstream ss;

    ss << dummy;

    ss >> alpha >> beta >> gamma >> delta >> value;


    int cant1 = cantor (alpha, beta);
    int cant2 = cantor (gamma, delta);

    int sign1 = sign_map.find(cant1)->second;
    int sign2 = sign_map.find(cant2)->second;

    int m1 = morten_map.find(cant1)->second;
    int m2 = morten_map.find(cant2)->second;

    h2_mat[m1][m2] = value * sign1 * sign2;

  }

  return;

}


/***************************************************************



***************************************************************/

int cantor (size_t a, size_t b)
{
	double mapped = 0.5 * (a + b) * (a + b + 1.0) + b;

	int value = (int) mapped;

	return value;
}