import os
import sys
import io
import tarfile
import libsbml
import glob
import tempfile
import logging
import shutil

sys.path.insert(0, '/home/')
import rpTool as rpReader
import rpToolCache

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)

'''
logging.disable(logging.INFO)
logging.disable(logging.WARNING)
'''

## RetroPath2.0 reader for local packages
#
#
def rp2Reader_mem(rpreader,
                  rp2_pathways,
                  rp2paths_pathways,
                  rp2paths_compounds,
                  upper_flux_bound,
                  lower_flux_bound,
                  maxRuleIds,
                  pathway_id,
                  compartment_id,
                  species_group_id,
                  pubchem_search,
                  outputTar):
    rpsbml_paths = rpreader.rp2ToSBML(rp2_pathways,
                                      rp2paths_pathways,
                                      rp2paths_compounds,
                                      None,
                                      upper_flux_bound,
                                      lower_flux_bound,
                                      maxRuleIds,
                                      pathway_id,
                                      compartment_id,
                                      pubchem_search,
                                      species_group_id)
    #pass the SBML results to a tar
    if rpsbml_paths=={}:
        logging.error('rpReader did not generate any results')
        return False
    #outputTar = io.BytesIO()
    #with open(outputTar, 'w:xz') as tf:
    with tarfile.open(fileobj=outputTar, mode='w:gz') as tf:
        for rpsbml_name in rpsbml_paths:
            data = libsbml.writeSBMLToString(rpsbml_paths[rpsbml_name].document).encode('utf-8')
            fiOut = io.BytesIO(data)
            info = tarfile.TarInfo(name=rpsbml_name)
            info.size = len(data)
            tf.addfile(tarinfo=info, fileobj=fiOut)
    return True


## RetroPath2.0 reader for local packages
#
#
def rp2Reader_hdd(rpreader,
                  rp2_pathways,
                  rp2paths_pathways,
                  rp2paths_compounds,
                  upper_flux_bound,
                  lower_flux_bound,
                  maxRuleIds,
                  pathway_id,
                  compartment_id,
                  species_group_id,
                  pubchem_search,
                  outputTar):
    # check that the files are not empty
    if sum(1 for line in open(rp2paths_compounds))<=1:
        logging.error('RP2paths compounds is empty')
        return False
    if sum(1 for line in open(rp2_pathways))<=1:
        logging.error('RP2 pathways is empty')
        return False
    if sum(1 for line in open(rp2paths_pathways))<=1:
        logging.error('RP2paths pathways is empty')
        return False
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        #Note the return here is {} and thus we can ignore it
        rpsbml_paths = rpreader.rp2ToSBML(rp2_pathways,
                                          rp2paths_pathways,
                                          rp2paths_compounds,
                                          tmpOutputFolder,
                                          upper_flux_bound,
                                          lower_flux_bound,
                                          maxRuleIds,
                                          pathway_id,
                                          compartment_id,
                                          species_group_id,
                                          pubchem_search)
        if len(glob.glob(tmpOutputFolder+'/*'))==0:
            logging.error('rpReader did not generate any results')
            return False
        with tarfile.open(outputTar, mode='w:gz') as ot:
            for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.rpsbml', '').replace('.xml', ''))+'.sbml.xml'
                info = tarfile.TarInfo(fileName)
                info.size = os.path.getsize(sbml_path)
                ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    return True

'''
##
#
#
def main_string(outputTar,
         upper_flux_bound=999999,
         lower_flux_bound=0,
         maxRuleIds=2,
         compartment_id='MNXC3',
         pathway_id='rp_pathway',
         species_group_id='central_species'):
        #pass the cache parameters to the rpReader
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
        #### MEM #####
        """
        if not rp2Reader_mem(rpreader,
                    rp2_pathways,
                    rp2paths_pathways,
                    rp2paths_compounds,
                    int(upper_flux_bound),
                    int(lower_flux_bound),
                    int(maxRuleIds),
                    pathway_id,
                    compartment_id,
                    species_group_id,
                    outputTar):
            abort(204)
        """
        #### HDD #####
        isOK = rp2Reader_hdd(rpreader,
                             rp2_pathways,
                             rp2paths_pathways,
                             rp2paths_compounds,
                             int(upper_flux_bound),
                             int(lower_flux_bound),
                             int(maxRuleIds),
                             pathway_id,
                             compartment_id,
                             species_group_id,
                             pubchem_search,
                             outputTar_bytes)
        if not isOK:
            logging.error('Function returned an error')
        ########IMPORTANT######
        outputTar_bytes.seek(0)
        #######################
        with open(outputTar, 'wb') as f:
            shutil.copyfileobj(outputTar_bytes, f, length=131072)
'''

##
#
#
def main_rp2(outputTar,
             rp2_pathways,
             rp2paths_pathways,
             rp2paths_compounds,
             upper_flux_bound=999999,
             lower_flux_bound=0,
             maxRuleIds=2,
             compartment_id='MNXC3',
             pathway_id='rp_pathway',
             species_group_id='central_species',
             pubchem_search=False):
        #pass the cache parameters to the rpReader
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
        rpreader.chebi_mnxm = rpcache.chebi_mnxm
        #outputTar_bytes = io.BytesIO()
        #### MEM #####
        """
        if not rp2Reader_mem(rpreader,
                    rp2_pathways,
                    rp2paths_pathways,
                    rp2paths_compounds,
                    int(upper_flux_bound),
                    int(lower_flux_bound),
                    int(maxRuleIds),
                    pathway_id,
                    compartment_id,
                    species_group_id,
                    outputTar):
            abort(204)
        """
        #### HDD #####
        isOK = rp2Reader_hdd(rpreader,
                             rp2_pathways,
                             rp2paths_pathways,
                             rp2paths_compounds,
                             int(upper_flux_bound),
                             int(lower_flux_bound),
                             int(maxRuleIds),
                             pathway_id,
                             compartment_id,
                             species_group_id,
                             pubchem_search,
                             outputTar)
        if not isOK:
            logging.error('Function returned an error')
        '''
        ########IMPORTANT######
        outputTar_bytes.seek(0)
        #######################
        with open(outputTar, 'wb') as f:
            shutil.copyfileobj(outputTar_bytes, f, length=131072)
        '''

##
#
#
def main_tsv(outputTar,
             tsvfile,
             upper_flux_bound=999999,
             lower_flux_bound=0,
             compartment_id='MNXC3',
             pathway_id='rp_pathway',
             species_group_id='central_species'):
        #pass the cache parameters to the rpReader
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
        rpreader.chebi_mnxm = rpcache.chebi_mnxm
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            rpreader.TSVtoSBML(tsvfile,
                               tmpOutputFolder,
                               upper_flux_bound,
                               lower_flux_bound,
                               compartment_id,
                               pathway_id,
                               species_group_id)
            logging.error(glob.glob(tmpOutputFolder+'/*'))
            logging.error(outputTar)
            if len(glob.glob(tmpOutputFolder+'/*'))==0:
                logging.error('rpReader did not generate any results')
                return False
            with tarfile.open(outputTar, mode='w:gz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.rpsbml', '').replace('.xml', ''))+'.sbml.xml'
                    info = tarfile.TarInfo(fileName)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
        return True
