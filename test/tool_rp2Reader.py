#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpReader REST service

"""
import argparse
import rpToolServe

##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to parse RP2 to generate rpSBML collection')
    parser.add_argument('-rp2paths_compounds', type=str)
    parser.add_argument('-rp2_pathways', type=str)
    parser.add_argument('-rp2paths_pathways', type=str)
    parser.add_argument('-maxRuleIds', type=str)
    parser.add_argument('-pathway_id', type=str)
    parser.add_argument('-compartment_id', type=str)
    parser.add_argument('-outputTar', type=str)
    params = parser.parse_args()
    rpToolServe.main(params.outputTar,
                     params.rp2paths_compounds,
                     params.rp2_pathways,
                     params.rp2paths_pathways,
                     params.maxRuleIds,
                     params.compartment_id,
                     params.pathway_id)
