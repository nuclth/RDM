#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 11:20:37 2017

@author: alex
"""


def create_tb_list (nmax, sp_list, tb_list):
    """Function to create the two-particle basis set. A double loop is run over 
    the sp orbital basis set taking unique pairs. 
    
    Arguments:
        nmax    - nmax truncation size
        sp_list - dataframe to hold single-particle orbital quantum numbers 
        tb_list - dataframe to hold two-particle basis states
        
    Returns:
        tb_list
            number - two particle number counter
            n1 - principal quantum number
            l1 - orbital quantum  number
            j1 - angular momentum quantum number 
            mj1 - ang. mom. quantum num. projection
            sp1 - sp orbital number
            and n2, l2, j2, mj2, sp2 for particle 2
    """
    
    # total number of sp states and counter
    total_orbs = len (sp_list)
    num = 1
    
    for orb1 in range (0, total_orbs):
        for orb2 in range (orb1+1, total_orbs):
            
            # sp values
            sp1 = sp_list.iloc[orb1].number            
            n1  = sp_list.iloc[orb1].n
            l1  = sp_list.iloc[orb1].l
            j1  = sp_list.iloc[orb1].j
            mj1 = sp_list.iloc[orb1].mj

            sp2 = sp_list.iloc[orb2].number
            n2  = sp_list.iloc[orb2].n
            l2  = sp_list.iloc[orb2].l
            j2  = sp_list.iloc[orb2].j
            mj2 = sp_list.iloc[orb2].mj

            # Ham. rotationally invariant so only accept states with total M = 0         
            if mj1 + mj2 != 0:
                continue
            
            ncount = 2*n1 + l1 + 2*n2 + l2
            
            if ncount > nmax:
                continue
            
            # convert from 2j and 2mj to j and mj
            mj1 *= 0.5
            mj2 *= 0.5
            
            j1 *= 0.5
            j2 *= 0.5
            
            # add basis state to dataframe
            tb_list = tb_list.append({'number': num, 
                                      'n1': n1, 'l1': l1, 'j1': j1, 'mj1': mj1,
                                      'n2': n2, 'l2': l2, 'j2': j2, 'mj2': mj2,
                                      'sp1' : sp1, 'sp2' : sp2
                                      }, ignore_index=True) 
    
            # increment counter
            num +=1
    
    return tb_list



def sort_tb (tb_list):
    """Function to sort the two-particle basis states by their parity. 
    Even parity states are given first and then odd. The final dataframe also
    adds a new column with the first two entries being the block sizes of
    the even and odd parity states (all other values in the column being 0).
    
    Argument:
        tb_list - dataframe to hold two-particle basis states
    """
    
    # create even/odd two-particle basis states
    even_parity = (tb_list['l1'] + tb_list['l2']) % 2 == 0
    odd_parity  = (tb_list['l1'] + tb_list['l2']) % 2 == 1

    # total number of states
    num_even = sum (a for a in even_parity)
    num_odd = sum (a for a in odd_parity)
    

    # create new dataframe for even/odd ordering
    shifted_list = tb_list[even_parity]
    shifted_list = shifted_list.append(tb_list[odd_parity])

    # reset index after moving states
    shifted_list = shifted_list.reset_index(drop=True)

    # add column for block sizes
    shifted_list['blocks'] = 0
    shifted_list['blocks'][0] = num_even
    shifted_list['blocks'][1] = num_odd


    return shifted_list
