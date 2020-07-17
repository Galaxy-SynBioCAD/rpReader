#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpReader REST service

"""
import sys
sys.path.insert(0, '/home/')
import argparse
import logging
import os
import rpToolServe


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to parse RP2 to generate rpSBML collection')
    parser.add_argument('-rp2_pathways', type=str)
    parser.add_argument('-rp2paths_pathways', type=str)
    parser.add_argument('-rp2paths_compounds', type=str)
    parser.add_argument('-output', type=str)
    parser.add_argument('-upper_flux_bound', type=int, default=999999)
    parser.add_argument('-lower_flux_bound', type=int, default=0)
    parser.add_argument('-maxRuleIds', type=int, default=2)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-compartment_id', type=str, default='MNXC3')
    parser.add_argument('-species_group_id', type=str, default='central_species')
    parser.add_argument('-sink_species_group_id', type=str, default='rp_sink_species')
    parser.add_argument('-pubchem_search', type=str, default='False')
    params = parser.parse_args()
    if params.maxRuleIds<=0:
        logging.error('Max rule ID cannot be less or equal than 0: '+str(params.maxRuleIds))
        exit(1)
    if params.pubchem_search=='True' or params.pubchem_search=='T' or params.pubchem_search=='true' or params.pubchem_search=='t':
        pubchem_search = True
    elif params.pubchem_search=='False' or params.pubchem_search=='F' or params.pubchem_search=='false' or params.pubchem_search=='f':
        pubchem_search = False
    else:
        logging.error('Cannot interpret pubchem_search input: '+str(params.pubchem_search))
        exit(1)
    if params.maxRuleIds<=0:
        logging.error('Max Rule ID cannot be less or equal than 0: '+str(params.maxRuleIds))
        exit(1)
    rpToolServe.main_rp2(params.output,
                         params.rp2_pathways,
                         params.rp2paths_pathways,
                         params.rp2paths_compounds,
                         int(params.upper_flux_bound),
                         int(params.lower_flux_bound),
                         int(params.maxRuleIds),
                         params.compartment_id,
                         params.pathway_id,
                         params.species_group_id,
                         params.sink_species_group_id,
                         pubchem_search)
