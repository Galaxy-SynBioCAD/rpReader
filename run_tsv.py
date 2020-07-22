#!/usr/bin/env python3
"""
Created on March 24 2020

@author: Melchior du Lac
@description: Convert tsv to SBML files

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
def main(tsvfile,
         output,
         upper_flux_bound=999999,
         lower_flux_bound=0,
         compartment_id='MNXC3',
         pathway_id='rp_pathway',
         species_group_id='central_species',
         sink_species_group_id='rp_sink_species'):
    docker_client = docker.from_env()
    image_str = 'brsynth/rpreader-standalone:v2'
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
        shutil.copy(tsvfile, tmpOutputFolder+'/tsvfile.tsv')
        command = ['python',
                   '/home/tool_tsvReader.py',
                   '-tsvfile',
                   '/home/tmp_output/tsvfile.tsv',
                   '-upper_flux_bound',
                   str(upper_flux_bound),
                   '-lower_flux_bound',
                   str(lower_flux_bound),
                   '-pathway_id',
                   str(pathway_id),
                   '-compartment_id',
                   str(compartment_id),
                   '-species_group_id',
                   str(species_group_id),
                   '-sink_species_group_id',
                   str(sink_species_group_id),
                   '-output',
                   '/home/tmp_output/output.dat']
        container = docker_client.containers.run(image_str,
                                                 command,
                                                 detach=True,
                                                 stderr=True,
                                                 volumes={tmpOutputFolder+'/': {'bind': '/home/tmp_output', 'mode': 'rw'}})
        container.wait()
        err = container.logs(stdout=False, stderr=True)
        err_str = err.decode('utf-8') 
        if not 'ERROR' in err_str:
            shutil.copy(tmpOutputFolder+'/output.dat', output)
        else:
            print(err_str)
        #shutil.copy(tmpOutputFolder+'/output.dat', output)
        container.remove()


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Convert TSV file to SBML files')
    parser.add_argument('-tsvfile', type=str)
    parser.add_argument('-upper_flux_bound', type=int, default=999999)
    parser.add_argument('-lower_flux_bound', type=int, default=0)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-compartment_id', type=str, default='MNXC3')
    parser.add_argument('-species_group_id', type=str, default='central_species')
    parser.add_argument('-sink_species_group_id', type=str, default='rp_sink_species')
    parser.add_argument('-output', type=str)
    params = parser.parse_args()
    main(params.tsvfile,
         params.output,
         params.upper_flux_bound,
         params.lower_flux_bound,
         params.compartment_id,
         params.pathway_id,
         params.species_group_id,
         params.sink_species_group_id)
