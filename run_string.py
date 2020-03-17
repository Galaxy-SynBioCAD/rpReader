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
        command = ['/home/single_string.py',
                   '-rp2path',
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
        docker_client.containers.run(image_str,
                command,
                auto_remove=True,
                detach=False,
                volumes={tmpOutputFolder+'/': {'bind': '/home/tmp_output', 'mode': 'rw'}})
        shutil.copy(tmpOutputFolder+'/output.dat', outputTar)


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Convert the results of RP2 and rp2paths to SBML files')
    parser.add_argument('-reacString', type=str)
    parser.add_argument('-ec', type=str)
    params = parser.parse_args()
    main(params.reacString,
         params.ec)
