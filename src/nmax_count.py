#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 11:03:52 2017

@author: alex

Program to create a list of single-particle (sp) orbitals, the associated
two-body (TB) basis set in m-scheme, and related flag files for use in the
SDP program. 

User input is to specify the nmax value (N truncation scheme) for the basis. The
program does run a bit slow if this starts to get large say around nmax = 6 or 8.
Some optimization required here for larger basis sets.

The program output is a number of text data files:
    sp_list - list of the different single-particle states
    no_terms - non-zero matrix entries of the OBDM used in the N and O flags
    tb_list - list of the unique two-body basis states
    h2_terms - non-zero matrix entries for the Hamiltonian and TBDM
    p_terms - non-zero matrix entries for the P flag
"""

# standard libraries import
import sys

# import third party libraries
import pandas as pd
import numpy as np

# import custom modules    
import sp_functs
import tb_functs
import flag_functs

if __name__ == '__main__':
        
    nmax = 2
    
    # create single-particle and two-particle dataframes
    sp_list = pd.DataFrame(columns= ['number', 'n', 'l', '2j', '2mj'])
    tb_list = pd.DataFrame(columns= ['number', 'n1', 'l1', 'j1', 'mj1', 
                                     'n2', 'l2', 'j2', 'mj2', 'sp1', 'sp2'])


    # populate sp dataframe
    sp_list = sp_functs.create_sp_list (nmax, sp_list)


    # recast values as ints
    sp_list['number'] = sp_list['number'].astype(int)
    sp_list['n'] = sp_list['n'].astype(int)
    sp_list['l'] = sp_list['l'].astype(int)
    sp_list['2j'] = sp_list['2j'].astype(int)
    sp_list['2mj'] = sp_list['2mj'].astype(int)


    # sort sp dataframe and add blocks column
    sp_list = sp_functs.sort_sp (sp_list)    
    sp_list = sp_functs.add_blocks (sp_list)
 

    # create sp relational dataframe and no flag dataframe
    sp_relational = sp_functs.sp_relational_db_morten (nmax, sp_list)
    no_terms = flag_functs.no_flag (sp_list)


    # save sp list and no flag dataframes to disk
    np.savetxt('../me_files/ref_files/nmax' + str(nmax) + 
               '_python_sp.dat', sp_list.values, fmt='%6i')
    np.savetxt('../flag_files/nmax' + str(nmax) + 
               '_python_noflag.dat', no_terms, fmt='%6i')


    print ("OB TERMS COMPLETE")
    

    
    # create and sort two particle basis    
    tb_list = tb_functs.create_tb_list (nmax, sp_list, tb_list)
    tb_list = tb_functs.sort_tb (tb_list)

    # recast values as ints
    tb_list.number = tb_list.number.astype(int)
    tb_list.n1 = tb_list.n1.astype(int)
    tb_list.l1 = tb_list.l1.astype(int)
    tb_list.n2 = tb_list.n2.astype(int)
    tb_list.l2 = tb_list.l2.astype(int)
    tb_list.sp1 = tb_list.sp1.astype(int)
    tb_list.sp2 = tb_list.sp2.astype(int)


    # save tb list dataframe to disk
    np.savetxt('../me_files/ref_files/nmax' + str(nmax) + 
               '_python_tb.dat', tb_list.values, fmt='%6.2f')
    
    print ("TB LIST COMPLETE")
    
   
#    print ('Memory size sp list {} bytes'.format(sys.getsizeof(sp_list)))
#    print ('Memory size tb list {} bytes'.format(sys.getsizeof(tb_list)))

    # create and sort flag for hamiltonian dataframe
    h2_terms = flag_functs.h2_flag (tb_list, sp_relational)
    h2_terms = h2_terms.astype(int)


    # create and sort p flag dataframe   
    p_terms = flag_functs.p_flag (no_terms, h2_terms)
    p_terms = flag_functs.sort_pflag (p_terms)


    # save h2 and p flag info to disk
    np.savetxt('../flag_files/nmax' + str (nmax) +
               '_python_h2flag.dat', h2_terms, fmt='%6i')
    np.savetxt('../flag_files/nmax' + str (nmax) +
               '_python_pflag.dat', p_terms, fmt='%6i')
    
    
    print ("TB FLAGS COMPLETE - PROGRAM FINISHED")
