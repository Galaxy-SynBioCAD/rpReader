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
def main(reaction_string,
          ec,
          upper_flux_bound,
          lower_flux_bound,
          pathway_id,
          compartment_id,
          species_group_id,
          output)
    docker_client = docker.from_env()
    image_str = 'brsynth/rpreader-standalone:dev'
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
        command = ['/home/single_string.py',
                   '-reaction_string',
                   reaction_string,
                   '-ec',
                   ec,
                   '-upper_flux_bound',
                   upper_flux_bound,
                   '-lower_flux_bound',
                   lower_flux_bound,
                   '-pathway_id',
                   pathway_id,
                   '-compartment_id',
                   compartment_id,
                   '-species_group_id',
                   species_group_id,
                   '-output',
                   output]
        container = docker_client.containers.run(image_str,
                                                 command,
                                                 detach=True,
                                                 stderr=True,
                                                 volumes={tmpOutputFolder+'/': {'bind': '/home/tmp_output', 'mode': 'rw'}})
        container.wait()
        err = container.logs(stdout=False, stderr=True)
        err_str = err.decode('utf-8')
        print(err_str)
        if not 'ERROR' in err_str:
            shutil.copy(tmpOutputFolder+'/output.dat', output)
        container.remove()


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Convert the results of RP2 and rp2paths to SBML files')
    parser.add_argument('-reaction_string', type=str)
    parser.add_argument('-ec', type=str, default='')
    parser.add_argument('-upper_flux_bound', type=int, default=999999)
    parser.add_argument('-lower_flux_bound', type=int, default=0)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-compartment_id', type=str, default='MNXC3')
    parser.add_argument('-species_group_id', type=str, default='central_species')
    parser.add_argument('-output', type=str)
    params = parser.parse_args()
    main(params.reaction_string,
          params.ec,
          params.upper_flux_bound,
          params.lower_flux_bound,
          params.pathway_id,
          params.compartment_id,
          params.species_group_id,
          params.output)
