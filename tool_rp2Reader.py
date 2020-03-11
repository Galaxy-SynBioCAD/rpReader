#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpReader REST service

"""
import sys
sys.path.insert(0, '/home/')
import argparse
import rpToolServe

##
#
#
'''
def main(outputTar,
         rp2paths_compounds,
         rp2_pathways,
         rp2paths_pathways,
         upper_flux_bound=999999,
         lower_flux_bound=0,
         maxRuleIds=2,
         compartment_id='MNXC3',
         pathway_id='rp_pathway',
         species_group_id='central_species'):

'''
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to parse RP2 to generate rpSBML collection')
    parser.add_argument('-rp2paths_compounds', type=str)
    parser.add_argument('-rp2_pathways', type=str)
    parser.add_argument('-rp2paths_pathways', type=str)
    parser.add_argument('-upper_flux_bound', type=int, default=999999)
    parser.add_argument('-lower_flux_bound', type=int, default=0)
    parser.add_argument('-maxRuleIds', type=int)
    parser.add_argument('-pathway_id', type=str)
    parser.add_argument('-compartment_id', type=str)
    parser.add_argument('-species_group_id', type=str)
    parser.add_argument('-outputTar', type=str)
    params = parser.parse_args()
    rpToolServe.main(params.outputTar,
                     params.rp2paths_compounds,
                     params.rp2_pathways,
                     params.rp2paths_pathways,
                     int(params.upper_flux_bound),
                     int(params.lower_flux_bound),
                     int(params.maxRuleIds),
                     params.compartment_id,
                     params.pathway_id,
                     params.species_group_id)
