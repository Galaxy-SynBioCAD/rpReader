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


## RetroPath2.0 reader for local packages
#
#
def rp2Reader_mem(rpreader,
                  rp2paths_compounds,
                  rp2_pathways,
                  rp2paths_pathways,
                  upper_flux_bound,
                  lower_flux_bound,
                  maxRuleIds,
                  pathway_id,
                  compartment_id,
                  species_group_id,
                  outputTar):
    rpsbml_paths = rpreader.rp2ToSBML(rp2paths_compounds,
                                      rp2_pathways,
                                      rp2paths_pathways,
                                      None,
                                      upper_flux_bound,
                                      lower_flux_bound,
                                      maxRuleIds,
                                      pathway_id,
                                      compartment_id,
                                      species_group_id)
    #pass the SBML results to a tar
    if rpsbml_paths=={}:
        return False
    #outputTar = io.BytesIO()
    #with open(outputTar, 'w:xz') as tf:
    with tarfile.open(fileobj=outputTar, mode='w:xz') as tf:
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
                  rp2paths_compounds,
                  rp2_pathways,
                  rp2paths_pathways,
                  upper_flux_bound,
                  lower_flux_bound,
                  maxRuleIds,
                  pathway_id,
                  compartment_id,
                  species_group_id,
                  outputTar):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        #Note the return here is {} and thus we can ignore it
        rpsbml_paths = rpreader.rp2ToSBML(rp2paths_compounds,
                                          rp2_pathways,
                                          rp2paths_pathways,
                                          tmpOutputFolder,
                                          upper_flux_bound,
                                          lower_flux_bound,
                                          maxRuleIds,
                                          pathway_id,
                                          compartment_id,
                                          species_group_id)
        if len(glob.glob(tmpOutputFolder+'/*'))==0:
            return False
        with tarfile.open(fileobj=outputTar, mode='w:xz') as ot:
            for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.rpsbml', '').replace('.xml', ''))+'.rpsbml.xml'
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
                    rp2paths_compounds,
                    rp2_pathways,
                    rp2paths_pathways,
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
                             rp2paths_compounds,
                             rp2_pathways,
                             rp2paths_pathways,
                             int(upper_flux_bound),
                             int(lower_flux_bound),
                             int(maxRuleIds),
                             pathway_id,
                             compartment_id,
                             species_group_id,
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
             rp2paths_compounds,
             rp2_pathways,
             rp2paths_pathways,
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
                    rp2paths_compounds,
                    rp2_pathways,
                    rp2paths_pathways,
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
                             rp2paths_compounds,
                             rp2_pathways,
                             rp2paths_pathways,
                             int(upper_flux_bound),
                             int(lower_flux_bound),
                             int(maxRuleIds),
                             pathway_id,
                             compartment_id,
                             species_group_id,
                             outputTar_bytes)
        if not isOK:
            logging.error('Function returned an error')
        ########IMPORTANT######
        outputTar_bytes.seek(0)
        #######################
        with open(outputTar, 'wb') as f:
            shutil.copyfileobj(outputTar_bytes, f, length=131072)

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
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            rpreader.TSVtoSBML(tsvfile,
                               tmpOutputFolder,
                               upper_flux_bound,
                               lower_flux_bound,
                               compartment,
                               pathway_id,
                               species_group_id)
            if len(glob.glob(tmpOutputFolder+'/*'))==0:
                return False
            with tarfile.open(outputTar, mode='w:xz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.rpsbml', '').replace('.xml', ''))+'.rpsbml.xml'
                    info = tarfile.TarInfo(fileName)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
        return True
