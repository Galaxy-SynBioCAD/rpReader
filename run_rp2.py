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
def main(rp2_pathways,
         rp2paths_pathways,
         rp2paths_compounds,
         output,
         rules_rall='None',
         compounds='None',
         upper_flux_bound=999999,
         lower_flux_bound=0,
         maxRuleIds=2,
         compartment_id='MNXC3',
         pathway_id='rp_pathway',
         species_group_id='central_species',
         sink_species_group_id='rp_sink_species',
         pubchem_search='False'):
    """Function parse the results of RetroPath2.0 and rp2paths including external rules file

    :param rp2_pathways: The RetroPath2.0 results scope file
    :param rp2paths_pathways: The rp2paths result pathway (out_paths) file
    :param rp2paths_compounds: The rp2paths result compounds file
    :param output: The output collection of rpSBML files
    :param rules_rall: The rules file (Default: None)
    :param compounds: The compound file (Default: None)
    :param upper_flux_bound: The default upper flux bound (Default: 999999)
    :param lower_flux_bound: The default lower flux bound (Default: 0)
    :param maxRuleIds: The maximal number of rules associated with each step (Default: 2)
    :param compartment_id: The compartment SBML id (Default: MNXC3)
    :param pathway_id: The Groups heterologous pathway id (Default: rp_pathway)
    :param species_group_id: The Groups id of the central species (Default: central_species)
    :param sink_species_group_id: The Groups id of the rp_sink_species (Default: rp_sink_species)
    :param pubchem_search: Use the pubchem database to search for missing cross reference (Default: False)

    :type rp2_pathways: str 
    :type rp2paths_pathways: str
    :type rp2paths_compounds: str
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
    :type pubchem_search: bool

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
        if os.path.exists(rp2paths_compounds) and os.path.exists(rp2paths_pathways) and os.path.exists(rp2_pathways):
            shutil.copy(rp2paths_compounds, tmpOutputFolder+'/rp2paths_compounds.csv')
            shutil.copy(rp2paths_pathways, tmpOutputFolder+'/rp2paths_pathways.csv')
            shutil.copy(rp2_pathways, tmpOutputFolder+'/rp2_pathways.csv')
            if os.path.exists(rules_rall) and os.path.exists(compounds):
                shutil.copy(rules_rall, tmpOutputFolder+'/rules_rall.tsv')
                shutil.copy(compounds, tmpOutputFolder+'/compounds.tsv')
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
                           '-rules_rall',
                           '/home/tmp_output/rules_rall.tsv',
                           '-compounds',
                           '/home/tmp_output/compounds.tsv',
                           '-maxRuleIds',
                           str(maxRuleIds),
                           '-pathway_id',
                           str(pathway_id),
                           '-compartment_id',
                           str(compartment_id),
                           '-species_group_id',
                           str(species_group_id),
                           '-sink_species_group_id',
                           str(sink_species_group_id),
                           '-pubchem_search',
                           str(pubchem_search),
                           '-output',
                           '/home/tmp_output/output.dat']
            else:
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
                           '-rules_rall',
                           'None',
                           '-compounds',
                           'None',
                           '-maxRuleIds',
                           str(maxRuleIds),
                           '-pathway_id',
                           str(pathway_id),
                           '-compartment_id',
                           str(compartment_id),
                           '-species_group_id',
                           str(species_group_id),
                           '-sink_species_group_id',
                           str(sink_species_group_id),
                           '-pubchem_search',
                           str(pubchem_search),
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
                shutil.copy(tmpOutputFolder+'/output.dat', output)
            container.remove()
        else:
            logging.error('Cannot find one or more of the input files: '+str(rp2paths_compounds))
            logging.error('Cannot find one or more of the input files: '+str(rp2paths_pathways))
            logging.error('Cannot find one or more of the input files: '+str(rp2_pathways))
            exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Convert the results of RP2 and rp2paths to SBML files')
    parser.add_argument('-rp2paths_compounds', type=str)
    parser.add_argument('-rp2_pathways', type=str)
    parser.add_argument('-rp2paths_pathways', type=str)
    parser.add_argument('-upper_flux_bound', type=int, default=999999)
    parser.add_argument('-lower_flux_bound', type=int, default=0)
    parser.add_argument('-rules_rall', type=str, default='None')
    parser.add_argument('-compounds', type=str, default='None')
    parser.add_argument('-maxRuleIds', type=int, default=2)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-compartment_id', type=str, default='MNXC3')
    parser.add_argument('-species_group_id', type=str, default='central_species')
    parser.add_argument('-sink_species_group_id', type=str, default='rp_sink_species')
    parser.add_argument('-pubchem_search', type=str, default='False')
    parser.add_argument('-output', type=str)
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
    main(params.rp2_pathways,
         params.rp2paths_pathways,
         params.rp2paths_compounds,
         params.output,
         params.rules_rall,
         params.compounds,
         params.upper_flux_bound,
         params.lower_flux_bound,
         params.maxRuleIds,
         params.compartment_id,
         params.pathway_id,
         params.species_group_id,
         params.sink_species_group_id,
         pubchem_search)
