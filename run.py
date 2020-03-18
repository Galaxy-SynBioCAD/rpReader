#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Extract the sink from an SBML into RP2 friendly format

"""
import argparse
import tempfile
import os
import logging
import shutil
import docker


##
#
#
def main(outputTar,
         rp2paths_compounds,
         rp2_pathways,
         rp2paths_pathways,
         upper_flux_bound,
         lower_flux_bound,
         maxRuleIds,
         compartment_id,
         pathway_id,
         species_group_id):
    docker_client = docker.from_env()
    image_str = 'brsynth/rpreader-standalone'
    try:
        image = docker_client.images.get(image_str)
    except docker.errors.ImageNotFound:
        logging.warning('Could not find the image, trying to pull it')
        try:
            docker_client.images.pull(image_str)
            image = docker_client.images.get(image_str)
        except docker.errors.ImageNotFound:
            logging.error('Cannot pull image: '+str(image_str))
            exit(1)
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        shutil.copy(rp2paths_compounds, tmpOutputFolder+'/rp2paths_compounds.csv')
        shutil.copy(rp2paths_pathways, tmpOutputFolder+'/rp2paths_pathways.csv')
        shutil.copy(rp2_pathways, tmpOutputFolder+'/rp2_pathways.csv')
        command = ['/home/tool_rp2Reader.py',
                   '-rp2paths_compounds',
                   '/home/tmp_output/rp2paths_compounds.csv',
                   '-rp2_pathways',
                   '/home/tmp_output/rp2_pathways.csv',
                   '-rp2paths_pathways',
                   '/home/tmp_output/rp2paths_pathways.csv',
                   '-upper_flux_bound',
                   str(upper_flux_bound),
                   '-lower_flux_bound',
                   str(lower_flux_bound),
                   '-maxRuleIds',
                   str(maxRuleIds),
                   '-pathway_id',
                   str(pathway_id),
                   '-compartment_id',
                   str(compartment_id),
                   '-species_group_id',
                   str(species_group_id),
                   '-outputTar',
                   '/home/tmp_output/output.dat']
        container = docker_client.containers.run(image_str,
												 command,
												 detach=True,
                                                 stderr=True,
												 volumes={tmpOutputFolder+'/': {'bind': '/home/tmp_output', 'mode': 'rw'}})
        container.wait()
        err = container.logs(stdout=False, stderr=True)
        print(err)
        shutil.copy(tmpOutputFolder+'/output.dat', outputTar)
        container.remove()


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Convert the results of RP2 and rp2paths to SBML files')
    parser.add_argument('-rp2paths_compounds', type=str)
    parser.add_argument('-rp2_pathways', type=str)
    parser.add_argument('-rp2paths_pathways', type=str)
    parser.add_argument('-upper_flux_bound', type=int, default=999999)
    parser.add_argument('-lower_flux_bound', type=int, default=0)
    parser.add_argument('-maxRuleIds', type=int, default=2)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-compartment_id', type=str, default='MNXC3')
    parser.add_argument('-species_group_id', type=str, default='central_species')
    parser.add_argument('-outputTar', type=str)
    params = parser.parse_args()
    if params.maxRuleIds<0:
        logging.error('Max Rule ID cannot be <0: '+str(params.maxRuleIds))
        exit(1)
    main(params.outputTar,
         params.rp2paths_compounds,
         params.rp2_pathways,
         params.rp2paths_pathways,
         params.upper_flux_bound,
         params.lower_flux_bound,
         params.maxRuleIds,
         params.compartment_id,
         params.pathway_id,
         params.species_group_id)
