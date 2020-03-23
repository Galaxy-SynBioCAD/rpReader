#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpReader REST service

"""
import sys
sys.path.insert(0, '/home/')
import argparse
import rpTool as rpReader
import rpToolCache
import logging

##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to parse RP2 to generate rpSBML collection')
    parser.add_argument('-reaction_string', type=str)
    parser.add_argument('-upper_flux_bound', type=int, default=999999)
    parser.add_argument('-lower_flux_bound', type=int, default=0)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-compartment_id', type=str, default='MNXC3')
    parser.add_argument('-species_group_id', type=str, default='central_species')
    parser.add_argument('-output', type=str)
    params = parser.parse_args()
    #make the tar.xz 
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        input_tar = tmpOutputFolder+'/tmp_input.tar.xz'
        output_tar = tmpOutputFolder+'/tmp_output.tar.xz'
        with tarfile.open(input_tar, mode='w:xz') as tf:
            #tf.add(params.input)
            info = tarfile.TarInfo('single.rpsbml.xml') #need to change the name since galaxy creates .dat files
            info.size = os.path.getsize(params.input)
            tf.addfile(tarinfo=info, fileobj=open(params.input, 'rb'))

        reacStr_to_sbml(self,
                        reac_sctring,
                        ec=[])



        rpreader = rpReader.rpReader()
        rpcache = rpToolCache.rpToolCache()
        rpreader.deprecatedMNXM_mnxm = rpcache.deprecatedMNXM_mnxm
        rpreader.deprecatedMNXR_mnxr = rpcache.deprecatedMNXR_mnxr
        rpreader.mnxm_strc = rpcache.mnxm_strc
        rpreader.inchikey_mnxm = rpcache.inchikey_mnxm
        rpreader.rr_reactions = rpcache.rr_reactions
        rpreader.chemXref = rpcache.chemXref
        rpreader.compXref = rpcache.compXref
        rpreader.nameCompXref = rpcache.nameCompXref
        outputTar_bytes = io.BytesIO()



        rpToolServe.main(input_tar,
                         params.full_sbml,
                         output_tar,
                         params.sim_type,
                         params.source_reaction,
                         params.target_reaction,
                         params.source_coefficient,
                         params.target_coefficient,
                         isMax,
                         params.fraction_of,
                         dontMerge,
                         params.pathway_id,
                         params.compartment_id)
        with tarfile.open(output_tar) as outTar:
            outTar.extractall(tmpOutputFolder)
        out_file = glob.glob(tmpOutputFolder+'/*.rpsbml.xml')
        if len(out_file)>1:
            logging.warning('There are more than one output file...')
        shutil.copy(out_file[0], params.output)




    if params.maxRuleIds<0:
        logging.error('Max rule ID cannot be less than 0: '+str(params.maxRuleIds))
        exit(1)
    rpToolServe.main_string(params.reaction_string,
                            int(params.upper_flux_bound),
                            int(params.lower_flux_bound),
                            params.compartment_id,
                            params.pathway_id,
                            params.species_group_id,
                            params.output)
