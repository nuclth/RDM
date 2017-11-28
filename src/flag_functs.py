#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 11:55:41 2017

@author: alex
"""

import pandas as pd

def no_flag (sp_list):
    """Function to create a dataframe that holds the non-zero matrix entries 
    for the one-body N and O flags. N flag ensures that OBDM traces to N and O
    flag ensures linear relations among OBDM. Dataframe columns include matrix indices, 
    sp numbers, and block number information for use in creating an SDP file.
    
    Argument: 
        sp_list - dataframe to hold sp orbital quantum numbers

    Returns:
        sp_vals - dataframe that holds the following:
            b1 - first matrix index for SDP output (starts at 1) 
            b2 - second matrix index for SDP output
            m1 - first conventional matrix index (starts at 0)
            m2 - second conventional matrix index
            sp1 - first sp orbital number 
            sp2 - second sp orbital number
            block - which matrix block the entry is in
    """
      
    sp_vals = pd.DataFrame(columns=['b1','b2','m1','m2','sp1','sp2', 'block'])

    block_list = pd.DataFrame(sp_list[['blocks']])
    
    # matrix index acess (m1, m2)
    m1 = -1
    m2 = 0 
  
    block = 0
  
    # loop over different OB blocks
    for num in block_list['blocks']:
        
        # skip sp entries until end of block
        if num < 1:
            continue

        # increment block counter
        block += 1


        # create loop over block's upper triangle row
        for val1 in range (1, num+1, 1):
            m1 += 1
            m2 = m1
            
            sp_num1 = sp_list.iloc[m1].number
            
            # loop over upper triangle column
            for val2 in range (val1, num+1, 1):
                sp_num2 = sp_list.iloc[m2].number
                sp_vals = sp_vals.append({'b1': val1, 'b2': val2, 
                                          'm1': m1, 'm2': m2, 'sp1' : sp_num1,
                                          'sp2' : sp_num2, 'block': block
                                          }, ignore_index=True)
                m2 += 1
                



    return sp_vals




def h2_flag (tb_list, sp_relational):
    """Function to create a dataframe that holds the non-zero matrix entries 
    for the h2 flag. The h2 flag gives the SDP style matrix index output for 
    the two-body part of the hamiltonian.
    
    Argument: 
        tb_list - dataframe to hold two-particle basis set
        sp_relational - dataframe that relates my sp numbers to Morten's

    Returns:
        sp_vals - dataframe that holds the following:
            b1 - first matrix index for SDP output (starts at 1) 
            b2 - second matrix index for SDP output
            m1 - first conventional matrix index (starts at 0)
            m2 - second conventional matrix index
            sp1-4 - my sp orbital numbers for the matrix element
            m_sp1-4 - Morten's sp orbital numbers for the matrix element
            block - which matrix block the entry is in
    """
    
    sp_vals = pd.DataFrame(columns=['b1','b2','m1','m2',
                                    'sp1','sp2','sp3','sp4', 
                                    'm_sp1', 'm_sp2', 'm_sp3', 'm_sp4',
                                    'block'])

    block_list = pd.DataFrame(tb_list[['blocks']])
    
    # matrix index access (m1, m2) 
    m1 = -1
    m2 = 0 
    block = 0
    
    
    # loop over two blocks: even and odd parity
    for num in block_list['blocks']:
        
        # skip null entries
        if num < 1:
            continue

        block += 1

        # create loop over block's upper triangle rows
        for val1 in range (1, num+1, 1):
            m1 += 1
            m2 = m1
            
            sp1 = tb_list.iloc[m1].sp1
            sp2 = tb_list.iloc[m1].sp2
            
            m_sp1 = sp_relational[sp_relational['alex_sp'] == sp1].morten_sp.iloc[0]
            m_sp2 = sp_relational[sp_relational['alex_sp'] == sp2].morten_sp.iloc[0]
            
            # loop over block's upper triangle columns
            for val2 in range (val1, num+1, 1):
                sp3 = tb_list.iloc[m2].sp1
                sp4 = tb_list.iloc[m2].sp2
                
                m_sp3 = sp_relational[sp_relational['alex_sp'] == sp3].morten_sp.iloc[0]
                m_sp4 = sp_relational[sp_relational['alex_sp'] == sp4].morten_sp.iloc[0]
                
                sp_vals = sp_vals.append({'b1': val1, 'b2': val2, 
                                          'm1': m1, 'm2': m2, 'sp1': sp1,
                                          'sp2': sp2, 'sp3':sp3, 
                                          'sp4' : sp4, 'm_sp1' : m_sp1, 'm_sp2': m_sp2,
                                          'm_sp3' : m_sp3, 'm_sp4': m_sp4, 'block': block
                                          }, ignore_index=True)
                m2 += 1
         

         
    return sp_vals


def p_flag (no_flag, h2_flag):
    """Function to create a dataframe that holds the non-zero matrix entries 
    for the p flag. P flag ensures linear relations between the OB and TB density
    matrices.
    
    Argument: 
        no_flag - dataframe to hold non-zero matrix entries for the OBDM
        h2_flag - dataframe to hold non-zero matrix entries for the TBDM

    Returns:
        sp_vals - dataframe that holds the following:
            ob_b1 - first matrix index for OB SDP output (starts at 1) 
            ob_b2 - second matrix index for OB SDP output
            ob_block - which OBDM block we're in
            tb_b1 - first matrix index for TB SDP output (starts at 1)
            tb_b2 - second matrix index for TB SDP output 
            tb_block - which TBDM block we're in (even or odd parity : 1 or 2)
            new - whether entry corresponds to new constraint matrix or not (True/False)
    """
    
    sp_vals = pd.DataFrame(columns=['ob_b1','ob_b2','ob_block',
                                    'tb_b1','tb_b2','tb_block','new'])
    
    # list to hold sp number repeats
    repeats = []
    
    # loop over first particle in TB matrix element
    for sp1, sp3 in zip(h2_flag.sp1, h2_flag.sp3):
        
        # skip if we've seen the sp numbers before
        if (sp1, sp3) in repeats:
            continue
        
        repeats.append((sp1,sp3))
        
       
        new = True
        
        # pull off OB row that matches sp numbers
        ob_row = no_flag[(no_flag['sp1'] == sp1) & (no_flag['sp2'] == sp3)]
        
        # check for row existence
        if ob_row.empty:
            continue
        
        # pull off SDP matrix indices and block info
        ob_b1 = ob_row.b1.iloc[0]
        ob_b2 = ob_row.b2.iloc[0]
        ob_block = ob_row.block.iloc[0]
        
        loop_repeats = []
        
        # loop over second particle in TB matrix element
        for sp2, sp4 in zip(h2_flag.sp2, h2_flag.sp4):
            
            # skip if trace condition not satisfied
            if (sp2 != sp4):
                continue

            # skip sp repeats
            if (sp2, sp4) in loop_repeats:
                continue
            
            
            loop_repeats.append((sp2,sp4))
            
            # pull off TB row that matches sp numbers
            tb_row = h2_flag[(h2_flag['sp1'] == sp1) & (h2_flag['sp2'] == sp2)
            & (h2_flag['sp3'] == sp3) & (h2_flag['sp4'] == sp4)]

            # check row existence
            if tb_row.empty:
                continue
            
            # pull off SDP matrix indices and block info
            tb_b1 = tb_row.b1.iloc[0]
            tb_b2 = tb_row.b2.iloc[0]
            tb_block = tb_row.block.iloc[0]


            sp_vals = sp_vals.append({'ob_b1':ob_b1, 'ob_b2':ob_b2, 
                                      'ob_block':ob_block, 'tb_b1':tb_b1, 'tb_b2':tb_b2,
                                      'tb_block':tb_block, 'new':new}, ignore_index=True)
            
            new = False
        

    return sp_vals


def sort_pflag (p_terms):
    """Function to sort our P flag terms by which constraint matrix we are in 
    (corresponding to the 'new' value), then each individual term by where we 
    are in the matrix ('tb_block' and 'tb_b1', 'tb_b2' values). 
    
    Argument: 
        p_terms - dataframe that holds non-zero matrix entries for P flag

    Returns:
        sp_vals - sorted p_terms dataframe (sorted by parity )
    """
    
    # final output dataframe
    sp_vals = pd.DataFrame(columns=['ob_b1','ob_b2','ob_block',
                                    'tb_b1','tb_b2','tb_block','new'])

    # intermediate dataframe holder
    hold_vals = pd.DataFrame(columns=['ob_b1','ob_b2','ob_block',
                                    'tb_b1','tb_b2','tb_block','new'])
    
    start = True   
    

    # loop over index and independent constraint matrix indicator
    for index, term in zip(p_terms.index, p_terms['new']):
        
        # append starting value and all entries in the same constraint matrix
        if term == False or start:
            hold_vals = hold_vals.append(p_terms.iloc[index], ignore_index=True)
        
        # sort terms in holder if we hit a new constraint matrix
        if term == True and not start:
            hold_vals = hold_vals.sort_values(['tb_block', 'tb_b1','tb_b2'])
            
            # reset sorted block 'new' values
            hold_vals.loc[:, 'new'] = False
            hold_vals.loc[0, 'new'] = True
            
            # pull off sorted constraint matrix block
            sp_vals = sp_vals.append(hold_vals, ignore_index=True)
            

            # drop terms from intermediate holder and append first entry of new matrix
            hold_vals.drop(hold_vals.index, inplace=True)
            hold_vals = hold_vals.append(p_terms.iloc[index], ignore_index=True)

        
        # sort for last entry in p_terms
        if index+1 == len(p_terms):
            hold_vals = hold_vals.sort_values(['tb_block', 'tb_b1','tb_b2'])
            hold_vals.loc[:, 'new'] = False
            hold_vals.loc[0, 'new'] = True
            sp_vals = sp_vals.append(hold_vals, ignore_index=True)
            
        start = False
     
    
    return sp_vals

