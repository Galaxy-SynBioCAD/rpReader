#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpReader REST service

"""
import requests
import argparse
import json


##
#
#
def rp2ReaderUpload(rp2paths_compounds,
        rp2_pathways,
        rp2paths_pathways,
        maxRuleIds,
        pathway_id,
        compartment_id,
        server_url,
        outputTar):
    # Post request
    data = {'maxRuleIds': maxRuleIds, 'pathway_id': pathway_id, 'compartment_id': compartment_id}
    files = {'rp2paths_compounds': open(rp2paths_compounds, 'rb'),
             'rp2paths_pathways': open(rp2paths_pathways, 'rb'),
             'rp2_pathways': open(rp2_pathways, 'rb'),
             'data': ('data.json', json.dumps(data))}
    r = requests.post(server_url+'/Query', files=files)
    r.raise_for_status()
    with open(outputTar, 'wb') as ot:
        ot.write(r.content)


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
    parser.add_argument('-server_url', type=str)
    parser.add_argument('-outputTar', type=str)
    params = parser.parse_args()
    rp2ReaderUpload(params.rp2paths_compounds,
            params.rp2_pathways,
            params.rp2paths_pathways,
            params.maxRuleIds,
            params.pathway_id,
            params.compartment_id,
            params.server_url,
            params.outputTar)
