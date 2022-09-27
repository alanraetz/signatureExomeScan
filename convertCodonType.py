# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 12:35:39 2022

@author: alanr
"""

Type = 'G[C>T]T'  # convert to GCT

codon = Type[0] + Type[2] + Type[6]


def get_wt_codon (string):
    return string[0] + string[2] + string[6]