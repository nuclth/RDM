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
        
        if num < 1:
            continue

        block += 1

        for val1 in range (1, num+1, 1):
            m1 += 1
            m2 = m1
            
            sp_num1 = sp_list.iloc[m1].number
            
            
            for val2 in range (val1, num+1, 1):
                sp_num2 = sp_list.iloc[m2].number
                sp_vals = sp_vals.append({'b1': val1, 'b2': val2, 
                                          'm1': m1, 'm2': m2, 'sp1' : sp_num1,
                                          'sp2' : sp_num2, 'block': block
                                          }, ignore_index=True)
                m2 += 1
                



    return sp_vals




def h2_flag (tb_list, sp_relational):
    
    sp_vals = pd.DataFrame(columns=['b1','b2','m1','m2',
                                    'sp1','sp2','sp3','sp4', 
                                    'm_sp1', 'm_sp2', 'm_sp3', 'm_sp4',
                                    'block'])

    block_list = pd.DataFrame(tb_list[['blocks']])
    
    m1 = -1
    m2 = 0 
    block = 0
    
    for num in block_list['blocks']:
        
        if num < 1:
            continue

        block += 1

        for val1 in range (1, num+1, 1):
            m1 += 1
            m2 = m1
            
            sp1 = tb_list.iloc[m1].sp1
            sp2 = tb_list.iloc[m1].sp2
            
            m_sp1 = sp_relational[sp_relational['alex_sp'] == sp1].morten_sp.iloc[0]
            m_sp2 = sp_relational[sp_relational['alex_sp'] == sp2].morten_sp.iloc[0]
            
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
         
#    print (sp_vals)
#    sys.exit()
         
    return sp_vals


def p_flag (no_flag, h2_flag):
    
    
#    sp_vals = pd.DataFrame({'sp1': tb_list.sp1, 'sp2':tb_list.sp2})

    sp_vals = pd.DataFrame(columns=['ob_b1','ob_b2','ob_block',
                                    'tb_b1','tb_b2','tb_block','new'])
    repeats = []
    

    for sp1, sp3 in zip(h2_flag.sp1, h2_flag.sp3):
        
        
        if (sp1, sp3) in repeats:
            continue
        
        repeats.append((sp1,sp3))
        
#        print (repeats)
        
        new = True

        

#        print ('SP1, SP3 {} {}'.format(sp1, sp3))
        
        ob_row = no_flag[(no_flag['sp1'] == sp1) & (no_flag['sp2'] == sp3)]
        
        if ob_row.empty:
            continue
        
#        print ('OB ROW')
#        print (ob_row)
        
        ob_b1 = ob_row.b1.iloc[0]
        ob_b2 = ob_row.b2.iloc[0]
        ob_block = ob_row.block.iloc[0]
        
        loop_repeats = []
        
        for sp2, sp4 in zip(h2_flag.sp2, h2_flag.sp4):
            
            if (sp2 != sp4):
                continue


            if (sp2, sp4) in loop_repeats:
                continue
            
            
            loop_repeats.append((sp2,sp4))
            
            tb_row = h2_flag[(h2_flag['sp1'] == sp1) & (h2_flag['sp2'] == sp2)
            & (h2_flag['sp3'] == sp3) & (h2_flag['sp4'] == sp4)]

            if tb_row.empty:
                continue
            
#            print ('TB_row')
#            print (tb_row)
            
            tb_b1 = tb_row.b1.iloc[0]
            tb_b2 = tb_row.b2.iloc[0]
            tb_block = tb_row.block.iloc[0]

#            print ('ENTRIES {} {} {} {} {} {} {}'.format(
#                    ob_b1, ob_b2, ob_block, tb_b1, tb_b2, tb_block, new))

            sp_vals = sp_vals.append({'ob_b1':ob_b1, 'ob_b2':ob_b2, 
                                      'ob_block':ob_block, 'tb_b1':tb_b1, 'tb_b2':tb_b2,
                                      'tb_block':tb_block, 'new':new}, ignore_index=True)
            
            new = False
        
#        print (no_flag)

        
#        for i2, sp3, sp4 in zip (sp_vals.index, sp_vals.sp1, sp_vals.sp2):
#            if sp2 == sp4 and i1 <= i2 and (sp1,sp3) not in repeats:
#                repeats.append((sp1, sp3))
#                print (sp1, sp2, sp3, sp4)
#                count+=1


    return sp_vals


def sort_pflag (p_terms):
    
    sp_vals = pd.DataFrame(columns=['ob_b1','ob_b2','ob_block',
                                    'tb_b1','tb_b2','tb_block','new'])

    hold_vals = pd.DataFrame(columns=['ob_b1','ob_b2','ob_block',
                                    'tb_b1','tb_b2','tb_block','new'])
    
    start = True   
    
#    p_terms['new'].iloc[len(p_terms)-1] = False

    for index, term in zip(p_terms.index, p_terms['new']):
        
        if term == False or start:
            hold_vals = hold_vals.append(p_terms.iloc[index], ignore_index=True)
        
#        print (index, term)
        
        if term == True and not start:
            hold_vals = hold_vals.sort_values(['tb_block', 'tb_b1','tb_b2'])
            hold_vals['new'] = False
            hold_vals['new'].iloc[0] = True
            sp_vals = sp_vals.append(hold_vals, ignore_index=True)
            hold_vals.drop(hold_vals.index, inplace=True)
            hold_vals = hold_vals.append(p_terms.iloc[index], ignore_index=True)

        
        if index+1 == len(p_terms):
            hold_vals = hold_vals.sort_values(['tb_block', 'tb_b1','tb_b2'])
            hold_vals['new'] = False
            hold_vals['new'].iloc[0] = True
            sp_vals = sp_vals.append(hold_vals, ignore_index=True)
            
        start = False
     
    
    return sp_vals

