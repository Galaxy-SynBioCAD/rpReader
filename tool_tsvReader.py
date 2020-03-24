#!/usr/bin/env python3
"""
Created on March 24 2019

@author: Melchior du Lac
@description: Galaxy script to parse TSV file to SBML

"""
import sys
sys.path.insert(0, '/home/')
import argparse
import rpToolServe
import logging

##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to parse RP2 to generate rpSBML collection')
    parser.add_argument('-tsvfile', type=str)
    parser.add_argument('-upper_flux_bound', type=int, default=999999)
    parser.add_argument('-lower_flux_bound', type=int, default=0)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-compartment_id', type=str, default='MNXC3')
    parser.add_argument('-species_group_id', type=str, default='central_species')
    parser.add_argument('-output', type=str)
    params = parser.parse_args()
    if params.maxRuleIds<0:
        logging.error('Max rule ID cannot be less than 0: '+str(params.maxRuleIds))
        exit(1)
    rpToolServe.main_tsv(params.output,
                         params.tsvfile,
                         int(params.upper_flux_bound),
                         int(params.lower_flux_bound),
                         params.compartment_id,
                         params.pathway_id,
                         params.species_group_id)

