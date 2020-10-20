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


def main(tsvfile,
         output,
         upper_flux_bound=999999,
         lower_flux_bound=0,
         compartment_id='MNXC3',
         pathway_id='rp_pathway',
         species_group_id='central_species',
         sink_species_group_id='rp_sink_species'):
    """Function parse the results of RetroPath2.0 and rp2paths including external rules file

    :param tsvfile: The path to the tsv file
    :param output: The output collection of rpSBML files
    :param upper_flux_bound: The default upper flux bound (Default: 999999)
    :param lower_flux_bound: The default lower flux bound (Default: 0)
    :param compartment_id: The compartment SBML id (Default: MNXC3)
    :param pathway_id: The Groups heterologous pathway id (Default: rp_pathway)
    :param species_group_id: The Groups id of the central species (Default: central_species)
    :param sink_species_group_id: The Groups id of the rp_sink_species (Default: rp_sink_species)

    :type tsvfile: str 
    :type output: str 
    :type rules_rall: str
    :type compounds: str
    :type upper_flux_bound: int
    :type lower_flux_bound: int
    :type maxRuleIds: int
    :type compartment_id: str
    :type pathway_id: str
    :type species_group_id: str
    :type sink_species_group_id: str

    :rtype: None
    :return: None
    """
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
        if os.path.exists(tsvfile):
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
            if 'ERROR' in err_str:
                print(err_str)
            elif 'WARNING' in err_str:
                print(err_str)
            if not os.path.exists(tmpOutputFolder+'/output.dat'):
                print('ERROR: Cannot find the output file: '+str(tmpOutputFolder+'/output.dat'))
            else:
                shutil.copy(tmpOutputFolder+'/output.dat', rp2paths_compounds)
            container.remove()
        else:
            logging.error('Cannot find one or more of the input files: '+str(tsvfile))
            exit(1)


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
