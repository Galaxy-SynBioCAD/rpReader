import csv
import os
import itertools
import pickle
import gzip
from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs
import sys
import random
import json
import copy
#from .setup_self.logger import self.logger
import logging
import io
import re
#import tarfile
import requests
import time

import rpSBML

class rpReader:
    """Call to transform the results of RetroPath2.0 and rp2paths to SBML files
    """
    def __init__(self):
        """Constructor of the class
        """
        self.logger = logging.getLogger(__name__)
        self.logger.debug('Starting instance of rpReader')
        self.deprecatedCID_cid = None
        self.deprecatedRID_rid = None
        self.cid_strc = None
        self.inchikey_cid = None
        self.rr_reactions = None
        self.cid_xref = None
        self.cid_name = None
        self.comp_xref = None
        self.xref_comp = None
        self.chebi_cid = None
        self.pubchem_inchi = {}
        self.pubchem_inchikey = {}
        self.pubchem_smiles = {}
        #####################
        #self.pubchem_sec_count = 0
        #self.pubchem_sec_start = 0.0
        self.pubchem_min_count = 0
        self.pubchem_min_start = 0.0


    #######################################################################
    ############################# PRIVATE FUNCTIONS #######################
    #######################################################################


    def _checkRIDdeprecated(self, rid):
        """Function to create a dictionnary of old to new reaction id's

        :param rid: Reaction ID

        :type rid: str

        :rtype: str
        :return: Reaction ID
        """
        try:
            return self.deprecatedRID_rid[rid]
        except KeyError:
            return rid


    def _checkCIDdeprecated(self, cid):
        """Function to create a dictionnary of old to new chemical id's

        :param rid: Compound ID

        :type rid: str

        :rtype: str
        :return: Compound ID
        """
        try:
            return self.deprecatedCID_cid[cid]
        except KeyError:
            return cid


    def _pubChemLimit(self):
        """Function to wait until the allowed number of requests can be made to pubChem

        No more than 5 requests per second.
        No more than 400 requests per minute.
        No longer than 300 second running time per minute.
        Requests exceeding limits are rejected (HTTP 503 error)

        :rtype: None
        :return: None
        """
        if self.pubchem_min_start==0.0:
            self.pubchem_min_start = time.time()
        #self.pubchem_sec_count += 1
        self.pubchem_min_count += 1
        #### requests per minute ####
        if self.pubchem_min_count>=500 and time.time()-self.pubchem_min_start<=60.0:
            self.logger.warning('Reached 500 requests per minute for pubchem... waiting a minute')
            time.sleep(60.0)
            self.pubchem_min_start = time.time()
            self.pubchem_min_count = 0
        elif time.time()-self.pubchem_min_start>60.0:
            self.pubchem_min_start = time.time()
            self.pubchem_min_count = 0


    def _pubchemStrctSearch(self, strct, itype='inchi'):
        """Try to retreive the xref from an inchi structure using pubchem

        :param strct: The input structure
        :param itype: The type of input. Valid options: inchi, inchikey, smiles

        :type strct: str
        :type itype: str

        :rtype: dict
        :return: The resulting cross reference and structures
        """
        self._pubChemLimit()
        try:
            r = requests.post('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/'+str(itype)+'/xrefs/SBURL/JSON', data={itype: strct})
            res_list = r.json()
        except json.decoder.JSONDecodeError:
            self.logger.warning('JSON decode error')
            return {}
        try:
            res_list = res_list['InformationList']['Information']
        except KeyError:
            self.logger.warning('pubchem JSON keyerror: '+str(res_list))
            return {}
        xref = {}
        if len(res_list)==1:
            self._pubChemLimit()
            try:
                prop = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+str(res_list[0]['CID'])+'/property/IUPACName,InChI,InChIKey,CanonicalSMILES/JSON')
                prop_list = prop.json()
            except json.decoder.JSONDecodeError:
                self.logger.warning('JSON decode error')
                return {}
            try:
                name = prop_list['PropertyTable']['Properties'][0]['IUPACName']
                inchi = prop_list['PropertyTable']['Properties'][0]['InChI']
                inchikey = prop_list['PropertyTable']['Properties'][0]['InChIKey']
                smiles = prop_list['PropertyTable']['Properties'][0]['CanonicalSMILES']
            except KeyError:
                self.logger.warning('pubchem JSON keyerror: '+str(prop_list))
                return {}
            #TODO: need to determine how long cobra cannot handle this
            #TODO: determine if names that are too long is the problem and if not remove this part
            if len(name)>30:
                self._pubChemLimit()
                try:
                    syn = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+str(res_list[0]['CID'])+'/synonyms/JSON')
                    syn_lst = syn.json()
                except json.decoder.JSONDecodeError:
                    self.logger.warning('pubchem JSON decode error')
                    return {}
                try:
                    syn_lst = syn_lst['InformationList']['Information'][0]['Synonym']
                    syn_lst = [x for x in syn_lst if not 'CHEBI' in x and not x.isupper()]
                    name = syn_lst[0] #need a better way instead of just the firs tone
                except KeyError:
                    self.logger.warning('pubchem JSON keyerror: '+str(syn.json()))
                    return {}
                except IndexError:
                    name = ''
            xref['pubchem'] = [str(res_list[0]['CID'])]
            for url in res_list[0]['SBURL']:
                if 'https://biocyc.org/compound?orgid=META&id=' in url:
                    if 'biocyc' not in xref:
                        xref['biocyc'] = []
                    xref['biocyc'].append(url.replace('https://biocyc.org/compound?orgid=META&id=', ''))
                if 'http://www.hmdb.ca/metabolites/' in url:
                    if 'hmdb' not in xref:
                        xref['hmdb'] = []
                    xref['hmdb'].append(url.replace('http://www.hmdb.ca/metabolites/', ''))
                if 'http://www.genome.jp/dbget-bin/www_bget?cpd:' in url:
                    if 'kegg_c' not in xref:
                        xref['kegg_c'] = []
                    xref['kegg_c'].append(url.replace('http://www.genome.jp/dbget-bin/www_bget?cpd:', ''))
                if 'http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:' in url:
                    if 'chebi' not in xref:
                        xref['chebi'] = []
                    xref['chebi'].append(url.replace('http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:', ''))
        elif len(res_list)==0:
            self.logger.warning('Could not find results for: '+str(strct))
            return {}
        else:
            self.logger.warning('There are more than one result for '+str(strct)+'... Ignoring')
            return {}
        return {'name': name, 'inchi': inchi, 'inchikey': inchikey, 'smiles': smiles, 'xref': xref}


    def _convert_depiction(self, idepic, itype='smiles', otype={'inchikey'}):
        """Convert chemical depiction to others type of depictions

        Usage example:
         - convert_depiction(idepic='CCO', otype={'inchi', 'smiles', 'inchikey'})
         - convert_depiction(idepic='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', itype='inchi', otype={'inchi', 'smiles', 'inchikey'})

        :param idepic: Input string
        :param itype: The type of input
        :param otype: Type of output. Valid options: inchi, smiles, inchikey

        :type idepic: str 
        :type itype: str
        :type otype: dict

        :rtype: dict
        :return: Dictionnary of results
        """
        # Import (if needed)
        if itype == 'smiles':
            rdmol = MolFromSmiles(idepic, sanitize=True)
        elif itype == 'inchi':
            rdmol = MolFromInchi(idepic, sanitize=True)
        else:
            raise NotImplementedError('"{}" is not a valid input type'.format(itype))
        if rdmol is None:  # Check imprt
            raise NotImplementedError('Import error from depiction "{}" of type "{}"'.format(idepic, itype))
        # Export
        odepic = dict()
        for item in otype:
            if item == 'smiles':
                odepic[item] = MolToSmiles(rdmol)  # MolToSmiles is tricky, one mays want to check the possible options..
            elif item == 'inchi':
                odepic[item] = MolToInchi(rdmol)
            elif item == 'inchikey':
                odepic[item] = MolToInchiKey(rdmol)
            else:
                raise NotImplementedError('"{}" is not a valid output type'.format(otype))
        return odepic


    ###############################################################
    ############################ RP2paths entry functions #########
    ###############################################################


    def rp2ToSBML(self,
                  rp2_pathways,
                  rp2paths_pathways,
                  rp2paths_compounds,
                  tmpOutputFolder=None,
                  upper_flux_bound=999999,
                  lower_flux_bound=0,
                  maxRuleIds=10,
                  pathway_id='rp_pathway',
                  compartment_id='MNXC3',
                  species_group_id='central_species',
                  sink_species_group_id='rp_sink_species',
                  pubchem_search=False):
        """Function to group all the functions for parsing RP2 output to SBML files
        
        Takes RP2paths's compounds.txt and out_paths.csv and RetroPaths's *_scope.csv files and generates SBML

        :param rp2_pathways: The RetroPath2.0 results scope file
        :param rp2paths_pathways: The rp2paths result pathway (out_paths) file
        :param rp2paths_compounds: The rp2paths result compounds file
        :param tmpOutputFolder: A folder to output the results (Default: None)
        :param upper_flux_bound: The default upper flux bound (Default: 999999)
        :param lower_flux_bound: The default lower flux bound (Default: 0)
        :param maxRuleIds: The maximal number of rules associated with each step (Default: 10)
        :param pathway_id: The Groups heterologous pathway id (Default: rp_pathway)
        :param compartment_id: The compartment SBML id (Default: MNXC3)
        :param species_group_id: The Groups id of the central species (Default: central_species)
        :param sink_species_group_id: The Groups id of the rp_sink_species (Default: rp_sink_species)
        :param pubchem_search: Use the pubchem database to search for missing cross reference (Default: False)

        :type rp2_pathways: str 
        :type rp2paths_pathways: str
        :type rp2paths_compounds: str
        :type tmpOutputFolder: str
        :type upper_flux_bound: int
        :type lower_flux_bound: int
        :type maxRuleIds: int
        :type pathway_id: str
        :type compartment_id: str
        :type species_group_id: str
        :type sink_species_group_id: str
        :type pubchem_search: bool

        :rtype: dict
        :return: Dictionnary of pathways results
        """
        self.logger.debug(maxRuleIds)
        rp_strc = self._compounds(rp2paths_compounds)
        rp_transformation, sink_molecules = self._transformation(rp2_pathways)
        return self._outPathsToSBML(rp_strc,
                                    rp_transformation,
                                    sink_molecules,
                                    rp2paths_pathways,
                                    upper_flux_bound,
                                    lower_flux_bound,
                                    tmpOutputFolder,
                                    maxRuleIds,
                                    pathway_id,
                                    compartment_id,
                                    species_group_id,
                                    sink_species_group_id,
                                    pubchem_search)


    def _compounds(self, path):
        """Function to parse the compounds.txt file
        
        Extract the smile and the structure of each compounds of RP2Path output. Method to parse all the RP output compounds.

        :param path: Path to the compounds file

        :type path: str 

        :rtype: dict
        :return: Dictionnary of compounds results
        """
        #self.rp_strc = {}
        rp_strc = {}
        try:
            #### we might pass binary in the REST version
            if isinstance(path, bytes):
                reader = csv.reader(io.StringIO(path.decode('utf-8')), delimiter='\t')
            else:
                reader = csv.reader(open(path, 'r', encoding='utf-8'), delimiter='\t')
            next(reader)
            for row in reader:
                self.logger.debug(row)
                rp_strc[row[0]] = {'smiles': row[1]}  #, 'structure':row[1].replace('[','').replace(']','')
                try:
                    rp_strc[row[0]]['inchi'] = self.cid_strc[row[0]]['inchi']
                except KeyError:
                    #try to generate them yourself by converting them directly
                    try:
                        resConv = self._convert_depiction(idepic=row[1], itype='smiles', otype={'inchi'})
                        rp_strc[row[0]]['inchi'] = resConv['inchi']
                    except NotImplementedError as e:
                        self.logger.warning('Could not convert the following SMILES to InChI: '+str(row[1]))
                try:
                    rp_strc[row[0]]['inchikey'] = self.cid_strc[row[0]]['inchikey']
                    #try to generate them yourself by converting them directly
                    #TODO: consider using the inchi writing instead of the SMILES notation to find the inchikey
                except KeyError:
                    try:
                        resConv = self._convert_depiction(idepic=row[1], itype='smiles', otype={'inchikey'})
                        rp_strc[row[0]]['inchikey'] = resConv['inchikey']
                    except NotImplementedError as e:
                        self.logger.warning('Could not convert the following SMILES to InChI key: '+str(row[1]))
        except (TypeError, FileNotFoundError) as e:
            self.logger.error('Could not read the compounds file ('+str(path)+')')
            raise RuntimeError
        return rp_strc


    def _transformation(self, path):
        """Function to parse the scope.csv file
        
        Extract the reaction rules from the retroPath2.0 output using the scope.csv file

        :param path: Path to the compounds file

        :type path: str 

        :rtype: tuple
        :return: The RetroPath transformation and the list of sink molecules
        """
        rp_transformation = {}
        sink_molecules = []
        #### we might pass binary in the REST version
        reader = None
        if isinstance(path, bytes):
            reader = csv.reader(io.StringIO(path.decode('utf-8')), delimiter=',')
        else:
            try:
                reader = csv.reader(open(path, 'r'), delimiter=',')
            except FileNotFoundError:
                self.logger.error('Could not read the compounds file: '+str(path))
                return {}
        next(reader)
        for row in reader:
            self.logger.debug(row)
            if not row[1] in rp_transformation:
                rp_transformation[row[1]] = {}
                rp_transformation[row[1]]['rule'] = row[2]
                rp_transformation[row[1]]['ec'] = [i.replace(' ', '') for i in row[11][1:-1].split(',') if not i.replace(' ', '')=='NOEC']
            if row[7]=='1':
                for i in row[8].replace(']', '').replace('[', '').replace(' ', '').split(','):
                    sink_molecules.append(i)
        self.logger.debug(rp_transformation)
        self.logger.debug(sink_molecules)
        return rp_transformation, list(set(sink_molecules))


    #TODO: make sure that you account for the fact that each reaction may have multiple associated reactions
    def _outPathsToSBML(self,
                        rp_strc,
                        rp_transformation,
                        sink_molecules,
                        rp2paths_outPath,
                        upper_flux_bound=999999,
                        lower_flux_bound=0,
                        tmpOutputFolder=None,
                        maxRuleIds=10,
                        pathway_id='rp_pathway',
                        compartment_id='MNXC3',
                        species_group_id='central_species',
                        sink_species_group_id='rp_sink_species',
                        pubchem_search=False):
        """Function to parse the out_paths.csv file

        Reading the RP2path output and extract all the information for each pathway RP2path Metabolic pathways from out_paths.csv create all the different values for heterologous paths from the RP2path out_paths.csv file. Note that path_step are in reverse order here

        :param rp_strc: The structures dictionary from RetroPath2.0
        :param rp_transformation: The transformation dictionary from RetroPath2.0
        :param sink_molecules: The dictionary of sink molecules from RetroPath2.0
        :param rp2paths_outPath: The rp2paths result file
        :param upper_flux_bound: The default upper flux bound (Default: 999999)
        :param lower_flux_bound: The default lower flux bound (Default: 0)
        :param tmpOutputFolder: A folder to output the results (Default: None)
        :param maxRuleIds: The maximal number of rules associated with each step (Default: 10)
        :param pathway_id: The Groups heterologous pathway id (Default: rp_pathway)
        :param compartment_id: The compartment SBML id (Default: MNXC3)
        :param species_group_id: The Groups id of the central species (Default: central_species)
        :param sink_species_group_id: The Groups id of the rp_sink_species (Default: rp_sink_species)
        :param pubchem_search: Use the pubchem database to search for missing cross reference (Default: False)

        :type rp2_strc: str 
        :type rp_transformation: str 
        :type sink_molecules: str 
        :type rp2paths_outPaths: str
        :type upper_flux_bound: int
        :type lower_flux_bound: int
        :type tmpOutputFolder: str
        :type maxRuleIds: int
        :type pathway_id: str
        :type compartment_id: str
        :type species_group_id: str
        :type sink_species_group_id: str
        :type pubchem_search: bool

        :rtype: dict
        :return: The dictionary of the rp2paths results or False boolean if fails
        """
        self.logger.debug(maxRuleIds)
        #try:
        rp_paths = {}
        sink_species = []
        #reactions = self.rr_reactionsingleRule.split('__')[1]s
        #with open(path, 'r') as f:
        #### we might pass binary in the REST version
        self.logger.debug('Parsing the following file: '+str(rp2paths_outPath))
        if isinstance(rp2paths_outPath, bytes):
            reader = csv.reader(io.StringIO(rp2paths_outPath.decode('utf-8')))
        else:
            reader = csv.reader(open(rp2paths_outPath, 'r'))
        next(reader)
        current_path_id = 0
        path_step = 1
        for row in reader:
            #in_sink
            #in_sink = int(row[7])
            self.logger.debug('Parsing the row: '+str(row))
            #Remove all illegal characters in SBML ids
            row[3] = row[3].replace("'", "").replace('-', '_').replace('+', '')
            try:
                if not int(row[0])==current_path_id:
                    path_step = 1
                else:
                    path_step += 1
                #important to leave them in order
                current_path_id = int(row[0])
            except ValueError:
                self.logger.error('Cannot convert path_id to int ('+str(row[0])+')')
                #return {}
                return False
            #################################
            ruleIds = row[2].split(',')
            if ruleIds==None:
                self.logger.warning('The rulesIds is None')
                #pass # or continue
                continue
            ###WARNING: This is the part where we select some rules over others
            # we do it by sorting the list according to their score and taking the topx
            tmp_rr_reactions = {}
            for r_id in ruleIds:
                try:
                    for rea_id in self.rr_reactions[r_id]:
                        tmp_rr_reactions[str(r_id)+'__'+str(rea_id)] = self.rr_reactions[r_id][rea_id]
                except KeyError:
                    self.logger.warning('Cannot find the following reaction rule: '+str(r_id)+'. Ignoring it...')
                    pass
            if len(ruleIds)>int(maxRuleIds):
                self.logger.warning('There are too many rules, limiting the number to top '+str(maxRuleIds))
                try:
                    ruleIds = [y for y,_ in sorted([(i, tmp_rr_reactions[i]['rule_score']) for i in tmp_rr_reactions])][:int(maxRuleIds)] 
                except KeyError:
                    self.logger.warning('Could not select topX due inconsistencies between rules ids and rr_reactions... selecting random instead')
                    ruleIds = random.sample(tmp_rr_reactions, int(maxRuleIds))
            else:
                ruleIds = tmp_rr_reactions
            sub_path_step = 1
            for singleRule in ruleIds:
                tmpReac = {'rule_id': singleRule.split('__')[0],
                           'rule_ori_reac': singleRule.split('__')[1],
                           'rule_score': self.rr_reactions[singleRule.split('__')[0]][singleRule.split('__')[1]]['rule_score'],
                           'right': {},
                           'left': {},
                           'path_id': int(row[0]),
                           'step': path_step,
                           'transformation_id': row[1][:-2]}
                ############ LEFT ##############
                for l in row[3].split(':'):
                    tmp_l = l.split('.')
                    try:
                        tmpReac['left'][self._checkCIDdeprecated(tmp_l[1])] = int(tmp_l[0])
                    except ValueError:
                        self.logger.error('Cannot convert tmp_l[0] to int ('+str(tmp_l[0])+')')
                        return False
                ############## RIGHT ###########
                for r in row[4].split(':'):
                    tmp_r = r.split('.')
                    try:
                        tmpReac['right'][self._checkCIDdeprecated(tmp_r[1])] = int(tmp_r[0])
                    except ValueError:
                        self.logger.error('Cannot convert tmp_r[0] to int ('+str(tmp_r[0])+')')
                        return False
                #################################
                if not int(row[0]) in rp_paths:
                    rp_paths[int(row[0])] = {}
                if not int(path_step) in rp_paths[int(row[0])]:
                    rp_paths[int(row[0])][int(path_step)] = {}
                rp_paths[int(row[0])][int(path_step)][int(sub_path_step)] = tmpReac
                sub_path_step += 1
        #### pathToSBML ####
        #------------- compartment --------------
        try:
            compid = self.xref_comp[compartment_id]
        except KeyError:
            self.logger.error('Could not Xref compartment_id ('+str(compartment_id)+')')
            return False
        #for pathNum in self.rp_paths:
        sbml_paths = {}
        for pathNum in rp_paths:
            #first level is the list of lists of sub_steps
            #second is itertools all possible combinations using product
            altPathNum = 1
            for comb_path in list(itertools.product(*[[(i,y) for y in rp_paths[pathNum][i]] for i in rp_paths[pathNum]])):
                steps = []
                for i, y in comb_path:
                    steps.append(rp_paths[pathNum][i][y])
                path_id = steps[0]['path_id']
                rpsbml = rpSBML.rpSBML('rp_'+str(path_id)+'_'+str(altPathNum))
                #1) create a generic Model, ie the structure and unit definitions that we will use the most
                ##### TODO: give the user more control over a generic model creation:
                #   -> special attention to the compartment
                rpsbml.genericModel('RetroPath_Pathway_'+str(path_id)+'_'+str(altPathNum),
                                    'RP_model_'+str(path_id)+'_'+str(altPathNum),
                                    self.comp_xref[compid],
                                    compartment_id,
                                    upper_flux_bound,
                                    lower_flux_bound)
                #2) create the pathway (groups)
                rpsbml.createPathway(pathway_id)
                rpsbml.createPathway(species_group_id)
                rpsbml.createPathway(sink_species_group_id)
                #3) find all the unique species and add them to the model
                all_meta = set([i for step in steps for lr in ['left', 'right'] for i in step[lr]])
                self.logger.debug('all_meta: '+str(all_meta))
                for meta in all_meta:
                    self.logger.debug('------- '+str(meta)+' --------')
                    spe_inchi = rp_strc[meta]['inchi']
                    spe_inchikey = rp_strc[meta]['inchikey']
                    spe_smiles = rp_strc[meta]['smiles']
                    chem_name = None
                    spe_xref = {}
                    try:
                        chem_name = self.cid_name[self._checkCIDdeprecated(meta)]
                    except KeyError:
                        #if you cannot find using cid, try to retreive it using its inchikey
                        try:
                            if 'inchikey' in rp_strc[meta]:
                                if rp_strc[meta]['inchikey']:
                                    self.logger.debug('Using InChIkey to find the name of the compound: '+str(rp_strc[meta]['inchikey']))
                                    # WARNING here we use MNX since as of now, there are only MNX data that is parsed correctly
                                    tmp_cids = [i for i in self.inchikey_cid[rp_strc[meta]['inchikey']] if i[:3]=='MNX']
                                    #TODO: handle multiple matches. For now we assume that multiple MNX means that there are deprecated versions of the tool
                                    if tmp_cids:
                                        chem_name = self.cid_name[self._checkCIDdeprecated(tmp_cids[0])]
                                        self.logger.debug('Found the name: '+str(chem_name))
                        except KeyError:
                            self.logger.warning('Cannot find the name for: '+str(meta))
                            chem_name = None
                    #compile as much info as you can
                    #xref
                    try:
                        spe_xref = self.cid_xref[self._checkCIDdeprecated(meta)]
                    except KeyError:
                        #if you cannot find using cid, try to retreive it using its inchikey
                        try:
                            if 'inchikey' in rp_strc[meta]:
                                if rp_strc[meta]['inchikey']:
                                    self.logger.info('Using the InChIkey to find the xref of this molecule: '+str(rp_strc[meta]['inchikey']))
                                    # WARNING here we use MNX since as of now, there are only MNX data that is parsed correctly
                                    tmp_cids = [i for i in self.inchikey_cid[rp_strc[meta]['inchikey']] if i[:3]=='MNX']
                                    self.logger.debug('tmp_cids: '+str(tmp_cids))
                                    #TODO: handle multiple matches. For now we assume that multiple MNX means that there are deprecated versions of the tool
                                    if tmp_cids:
                                        spe_xref = self.cid_xref[self._checkCIDdeprecated(tmp_cids[0])]
                                        self.logger.debug('Found the xref: '+str(spe_xref))
                        except KeyError:
                            self.logger.warning('Cannot find the xref for: '+str(meta))
                            spe_xref = {}
                    ###### Try to recover the structures ####
                    pubchem_inchi = None
                    pubchem_inchikey = None
                    pubchem_smiles = None
                    pubchem_xref = {}
                    #inchi
                    try:
                        if not spe_xref and pubchem_search and spe_inchi:
                            try:
                                pubchem_inchi = self.pubchem_inchi[spe_inchi]['inchi']
                                pubchem_inchikey = self.pubchem_inchi[spe_inchi]['inchikey']
                                pubchem_smiles = self.pubchem_inchi[spe_inchi]['smiles']
                                pubchem_xref = self.pubchem_inchi[spe_inchi]['xref'] 
                            except KeyError:
                                pubres = self._pubchemStrctSearch(spe_inchi, 'inchi')
                                if not chem_name:
                                    chem_name = pubres['name']
                                if 'chebi' in pubres['xref']:
                                    try:
                                        spe_xref = self.cid_xref[self.chebi_cid[pubres['xref']['chebi'][0]]]
                                    except KeyError:
                                        self.logger.warning('Failed to retreive the pubchem info using the InChiKey') 
                                if not pubchem_xref:
                                    pubchem_xref = pubres['xref']
                                if not pubchem_inchikey:
                                    pubchem_inchikey = pubres['inchikey']
                                if not pubchem_smiles:
                                    pubchem_smiles = pubres['smiles']
                    except KeyError:
                        self.logger.warning('Bad results from pubchem results')
                        self.pubchem_inchi[spe_inchi] = {}
                        pass
                    #inchikey
                    try:
                        if not spe_xref and pubchem_search and spe_inchikey:
                            try:
                                pubchem_inchi = self.pubchem_inchikey[spe_inchikey]['inchi']
                                pubchem_inchikey = self.pubchem_inchikey[spe_inchikey]['inchikey']
                                pubchem_smiles = self.pubchem_inchikey[spe_inchikey]['smiles']
                                pubchem_xref = self.pubchem_inchikey[spe_inchikey]['xref'] 
                            except KeyError:
                                #if not self.pubchem_inchikey[spe_inchikey]=={}:
                                pubres = self._pubchemStrctSearch(spe_inchikey, 'inchikey')
                                if not chem_name:
                                    chem_name = pubres['name']
                                if 'chebi' in pubres['xref']:
                                    try:
                                        spe_xref = self.cid_xref[self.chebi_cid[pubres['xref']['chebi'][0]]]
                                    except KeyError:
                                        self.logger.warning('Failed to retreive the pubchem info using the InChiKey') 
                                if not pubchem_xref:
                                    pubchem_xref = pubres['xref']
                                if not pubchem_inchi:
                                    pubchem_inchi = pubres['inchi']
                                if not pubchem_smiles:
                                    pubchem_smiles = pubres['smiles']
                    except KeyError:
                        self.logger.warning('Bad results from pubchem results')
                        self.pubchem_inchikey[spe_inchikey] = {}
                        pass
                    #smiles
                    try:
                        if not spe_xref and pubchem_search and spe_smiles:
                            try:
                                pubchem_inchi = self.pubchem_smiles[spe_smiles]['inchi']
                                pubchem_inchikey = self.pubchem_smiles[spe_smiles]['inchikey']
                                pubchem_smiles = self.pubchem_smiles[spe_smiles]['smiles']
                                pubchem_xref = self.pubchem_smiles[spe_smiles]['xref'] 
                            except KeyError:
                                #if not self.pubchem_smiles[spe_smiles]=={}:
                                pubres = self._pubchemStrctSearch(spe_smiles, 'smiles')
                                if not chem_name:
                                    chem_name = pubres['name']
                                if 'chebi' in pubres['xref']:
                                    try:
                                        spe_xref = self.cid_xref[self.chebi_cid[pubres['xref']['chebi'][0]]]
                                    except KeyError:
                                        self.logger.warning('Failed to retreive the pubchem info using the SMILES')
                                if not pubchem_xref:
                                    pubchem_xref = pubres['xref']
                                if not pubchem_inchi:
                                    pubchem_inchi = pubres['inchi']
                                if not pubchem_inchikey:
                                    pubchem_inchikey = pubres['inchikey']
                    except KeyError:
                        self.pubchem_smiles[spe_smiles] = {}
                        self.logger.warning('Bad results from pubchem results')
                        pass
                    if not spe_inchi:
                        spe_inchi = pubchem_inchi
                    if not spe_inchikey:
                        spe_inchikey = pubchem_inchikey
                    if not spe_smiles:
                        spe_smiles = pubchem_smiles
                    if pubchem_inchi:
                        self.pubchem_inchi[pubchem_inchi] = {'inchi': pubchem_inchi, 'smiles': pubchem_smiles, 'inchikey': pubchem_inchikey, 'xref': pubchem_xref}
                    if pubchem_inchikey:
                        self.pubchem_inchikey[pubchem_inchikey] = {'inchi': pubchem_inchi, 'smiles': pubchem_smiles, 'inchikey': pubchem_inchikey, 'xref': pubchem_xref}
                    if pubchem_smiles:
                        self.pubchem_smiles[pubchem_smiles] = {'inchi': pubchem_inchi, 'smiles': pubchem_smiles, 'inchikey': pubchem_inchikey, 'xref': pubchem_xref}
                    if not spe_xref:
                        spe_xref = pubchem_xref
                    #pass the information to create the species
                    if chem_name:
                        chem_name = chem_name.replace("'", "")
                    self.logger.debug('Creating species: '+str(chem_name)+' ('+str(meta)+')')
                    if meta in sink_molecules:
                        self.logger.debug('Species is sink: '+str(sink_species_group_id))
                        rpsbml.createSpecies(meta,
                                             compartment_id,
                                             chem_name,
                                             spe_xref,
                                             spe_inchi,
                                             spe_inchikey,
                                             spe_smiles,
                                             species_group_id,
                                             sink_species_group_id)
                    else:
                        rpsbml.createSpecies(meta,
                                             compartment_id,
                                             chem_name,
                                             spe_xref,
                                             spe_inchi,
                                             spe_inchikey,
                                             spe_smiles,
                                             species_group_id)
                #4) add the complete reactions and their annotations
                for step in steps:
                    self.logger.debug('#########################################')
                    self.logger.debug(step)
                    self.logger.debug('#########################################')
                    #add the substep to the model
                    step['sub_step'] = altPathNum
                    self.logger.debug('Creating reaction: '+str('RP'+str(step['step'])))
                    self.logger.debug('Steps:'+str(step))
                    rpsbml.createReaction('RP'+str(step['step']), # parameter 'name' of the reaction deleted : 'RetroPath_Reaction_'+str(step['step']),
                                          upper_flux_bound,
                                          lower_flux_bound,
                                          step,
                                          compartment_id,
                                          rp_transformation[step['transformation_id']]['rule'],
                                          {'ec': rp_transformation[step['transformation_id']]['ec']},
                                          pathway_id)
                #5) adding the consumption of the target
                targetStep = {'rule_id': None,
                              'left': {[i for i in all_meta if i[:6]=='TARGET'][0]: 1},
                              'right': [],
                              'step': None,
                              'sub_step': None,
                              'path_id': None,
                              'transformation_id': None,
                              'rule_score': None,
                              'rule_ori_reac': None}
                self.logger.debug('Creating reaction: RP1_sink')
                rpsbml.createReaction('RP1_sink',
                                      upper_flux_bound,
                                      lower_flux_bound,
                                      targetStep,
                                      compartment_id)
                #6) Add the flux objectives
                if tmpOutputFolder:
                    rpsbml.writeSBML(tmpOutputFolder)
                else:
                    sbml_paths['rp_'+str(step['path_id'])+'_'+str(altPathNum)] = rpsbml
                altPathNum += 1
        return sbml_paths



    #############################################################################################
    ############################### TSV data tsv ################################################
    #############################################################################################


    #TODO: update this to the new compartements and others
    def _parseTSV(self, inFile, remove_inchi_4p=False, mnxHeader=False):
        """Function to parse the TSV of measured heterologous pathways to SBML

        Given the TSV of measured pathways, parse them to a dictionnary, readable to next be parsed to SBML

        :param inFile: The input TSV file
        :param remove_inchi_4p: Remove the 4th dimentsion of inchi's when adding them to SBML files
        :param mnxHeader: Reorganise the results around the target MNX products

        :type inFile: str
        :type remove_inchi_4p: bool
        :type mnxHeader: bool

        :rtype: dict
        :return: The dictionary of the pathways
        """
        data = {}
        try:
            for row in csv.DictReader(open(inFile), delimiter='\t'):
                ######## path_id ######
                try:
                    pathID = int(row['pathway_ID'])
                except ValueError:
                    self.logger.error('Cannot convert pathway ID: '+str(row['pathway_ID']))
                    continue
                if not pathID in data:
                    data[pathID] = {}
                    data[pathID]['isValid'] = True
                    data[pathID]['steps'] = {}
                ####### target #########
                if not 'target' in data[pathID]:
                    data[pathID]['target'] = {}
                    data[pathID]['target']['name'] = row['target_name']
                    if remove_inchi_4p:
                        data[pathID]['target']['inchi'] = '/'.join([row['target_structure'].split('/')[i] for i in range(len(row['target_structure'].split('/'))) if i<4])
                    else:
                        data[pathID]['target']['inchi'] = row['target_structure']
                ####### step #########
                try:
                    stepID = int(row['step'])
                except ValueError:
                    self.logger.error('Cannot convert step ID: '+str(row['step']))
                    data[pathID]['isValid'] = False
                    continue
                if stepID==0:
                    continue
                elif stepID==1:
                    data[pathID]['organism'] = row['organism'].replace(' ', '')
                    data[pathID]['reference'] = row['reference'].replace(' ', '')
                data[pathID]['steps'][stepID] = {}
                ##### substrates #########
                data[pathID]['steps'][stepID]['substrates'] = []
                lenDBref = len(row['substrate_dbref'].split(';'))
                for i in row['substrate_dbref'].split(';'):
                    if i=='':
                        lenDBref -= 1
                lenStrc = len(row['substrate_structure'].split('_'))
                for i in row['substrate_structure'].split('_'):
                    if i=='':
                        lenStrc -= 1
                lenSub = len(row['substrate_name'].split(';'))
                for i in row['substrate_name'].split(';'):
                    if i=='':
                        lenSub -= 1
                if lenSub==lenStrc==lenSub:
                    for name, inchi, dbrefs in zip(row['substrate_name'].split(';'),
                            row['substrate_structure'].split('_'),
                            row['substrate_dbref'].split(';')):
                        tmp = {}
                        if remove_inchi_4p:
                            tmp['inchi'] = '/'.join([inchi.split('/')[i] for i in range(len(inchi.split('/'))) if i<4])
                        else:
                            tmp['inchi'] = inchi.replace(' ', '')
                        tmp['name'] = name
                        tmp['dbref'] = {}
                        for dbref in dbrefs.split('|'):
                            if len(dbref.split(':'))==2:
                                db_name = dbref.split(':')[0].replace(' ', '').lower()
                                db_cid = dbref.split(':')[1].replace(' ', '')
                                if not db_name in tmp['dbref']:
                                    tmp['dbref'][db_name] = []
                                tmp['dbref'][db_name].append(db_cid)
                            else:
                                self.logger.warning('Ignoring the folowing product dbref ('+str(name)+'): '+str(dbref))
                                data[pathID]['isValid'] = False
                        data[pathID]['steps'][stepID]['substrates'].append(tmp)
                else:
                    self.logger.warning('Not equal length between substrate names, their structure or dbref ('+str(name)+'): '+str(row['substrate_name'])+' <--> '+str(row['substrate_structure'])+' <--> '+str(row['substrate_dbref']))
                    data[pathID]['isValid'] = False
                    continue
                ##### products #########
                data[pathID]['steps'][stepID]['products'] = []
                lenDBref = len(row['product_dbref'].split(';'))
                for i in row['product_dbref'].split(';'):
                    if i=='':
                        lenDBref -= 1
                lenStrc = len(row['product_structure'].split('_'))
                for i in row['product_structure'].split('_'):
                    if i=='':
                        lenStrc -= 1
                lenSub = len(row['product_name'].split(';'))
                for i in row['product_name'].split(';'):
                    if i=='':
                        lenSub -= 1
                if lenSub==lenStrc==lenDBref:
                    for name, inchi, dbrefs in zip(row['product_name'].split(';'),
                            row['product_structure'].split('_'),
                            row['product_dbref'].split(';')):
                        tmp = {}
                        if remove_inchi_4p:
                            tmp['inchi'] = '/'.join([inchi.split('/')[i] for i in range(len(inchi.split('/'))) if i<4])
                        else:
                            tmp['inchi'] = inchi.replace(' ', '')
                        tmp['name'] = name
                        tmp['dbref'] = {}
                        for dbref in dbrefs.split('|'):
                            if len(dbref.split(':'))==2:
                                db_name = dbref.split(':')[0].replace(' ', '').lower()
                                db_cid = dbref.split(':')[1].replace(' ', '')
                                if not db_name in tmp['dbref']:
                                    tmp['dbref'][db_name] = []
                                tmp['dbref'][db_name].append(db_cid)
                            else:
                                data[pathID]['isValid'] = False
                                self.logger.warning('Ignoring the folowing product dbref ('+str(name)+'): '+str(dbref))
                        data[pathID]['steps'][stepID]['products'].append(tmp)
                else:
                    self.logger.warning('Not equal length between substrate names, their structure or dbref ('+str(name)+'): '+str(row['product_name'])+' <--> '+str(row['product_structure'])+' <--> '+str(row['product_dbref']))
                    data[pathID]['isValid'] = False
                if not row['uniprot']=='':
                    data[pathID]['steps'][stepID]['uniprot'] = row['uniprot'].replace(' ', '').split(';')
                if not row['EC_number']=='':
                    data[pathID]['steps'][stepID]['ec_numbers'] = [i.replace(' ', '') for i in row['EC_number'].split(';')]
                data[pathID]['steps'][stepID]['enzyme_id'] = [i.replace(' ', '') for i in row['enzyme_identifier'].split(';')]
                data[pathID]['steps'][stepID]['enzyme_name'] = row['enzyme_name'].split(';')
        except FileNotFoundError:
            self.logger.error('Cannot open the file: '+str(inFile))
        #now loop through all of them and remove the invalid paths
        toRet = copy.deepcopy(data)
        for path_id in data.keys():
            if toRet[path_id]['isValid']==False:
                del toRet[path_id]
            else:
                del toRet[path_id]['isValid']
        #reorganise the results around the target products mnx
        if not mnxHeader:
            return toRet
        else:
            toRetTwo = {}
            for path_id in toRet:
                try:
                    final_pro_mnx = toRet[path_id]['steps'][max(toRet[path_id]['steps'])]['products'][0]['dbref']['mnx'][0]
                except KeyError:
                    self.logger.error('The species '+str(toRet[path_id]['steps'][max(toRet[path_id]['steps'])]['products'][0]['name'])+' does not contain a mnx database reference... skipping whole pathway number '+str(path_id))
                    #continue
                if not final_pro_mnx in toRetTwo:
                    toRetTwo[final_pro_mnx] = {}
                toRetTwo[final_pro_mnx][path_id] = toRet[path_id]
            return toRetTwo


    # TODO: update this with the new SBML groups (sink species) -- perhaps not applicable
    def TSVtoSBML(self,
                  inFile,
                  tmpOutputFolder=None,
                  upper_flux_bound=99999,
                  lower_flux_bound=0,
                  compartment_id='MNXC3',
                  pathway_id='rp_pathway',
                  species_group_id='central_species',
                  header_name=''):
        """Parse the TSV file to SBML format and adds them to the self.sbml_paths

        :param inFile: The input TSV file
        :param tmpOutputFolder: A folder to output the results (Default: None)
        :param upper_flux_bound: The default upper flux bound (Default: 999999)
        :param lower_flux_bound: The default lower flux bound (Default: 0)
        :param compartment_id: The compartment SBML id (Default: MNXC3)
        :param pathway_id: The Groups heterologous pathway id (Default: rp_pathway)
        :param species_group_id: The Groups id of the central species (Default: central_species)
        :param header_name: Overwrite the name given to the SBML files generated

        :type inFile: str
        :type tmpOutputFolder: str
        :type upper_flux_bound: int
        :type lower_flux_bound: int
        :type compartment_id: str
        :type pathway_id: str
        :type species_group_id: str
        :type header_name: str

        :rtype: dict
        :return: The dictionary of the pathways
        """
        data = self._parseTSV(inFile)
        sbml_paths = {}
        if header_name=='':
            header_name = inFile.split('/')[-1].replace('.tsv', '').replace('.csv', '')
        #TODO: need to exit at this loop
        for path_id in data:
            try:
                mnxc = self.xref_comp[compartment_id]
            except KeyError:
                self.logger.error('Could not Xref compartment_id ('+str(compartment_id)+')')
                return False
            rpsbml = rpSBML.rpSBML(header_name+'_'+str(path_id))
            #1) create a generic Model, ie the structure and unit definitions that we will use the most
            ##### TODO: give the user more control over a generic model creation:
            #   -> special attention to the compartment
            rpsbml.genericModel(header_name+'_Path'+str(path_id),
                                header_name+'_Path'+str(path_id),
                                self.comp_xref[mnxc],
                                compartment_id,
                                upper_flux_bound,
                                lower_flux_bound)
            #2) create the pathway (groups)
            rpsbml.createPathway(pathway_id)
            rpsbml.createPathway(species_group_id)
            #3) find all the unique species and add them to the model
            allChem = []
            for stepNum in data[path_id]['steps']:
                #because of the nature of the input we need to remove duplicates
                for i in data[path_id]['steps'][stepNum]['substrates']+data[path_id]['steps'][stepNum]['products']:
                    if not i in allChem:
                        allChem.append(i)            
            #add them to the SBML
            for chem in allChem:
                #PROBLEM: as it stands one expects the meta to be MNX
                if 'mnx' in chem['dbref']:
                    #must list the different models
                    meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                else:
                    #TODO: add the species with other types of xref in annotation
                    self.logger.warning('Some species are not referenced by a MNX id and will be ignored')
                    #try CHEBI
                    try:
                        meta = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                        meta = 'CHEBI_'+str(meta)
                    except KeyError:
                        #TODO: need to find a better way
                        self.logger.warning('Cannot determine MNX or CHEBI entry, using random')
                        tmpDB_name = list(chem['dbref'].keys())[0]
                        meta = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                        meta = str(tmpDB_name)+'_'+str(meta) 
                    #break
                #try to conver the inchi into the other structures
                smiles = None
                inchikey = None
                try:
                    resConv = self._convert_depiction(idepic=chem['inchi'], itype='inchi', otype={'smiles','inchikey'})
                    smiles = resConv['smiles']
                    inchikey = resConv['inchikey']
                except NotImplementedError as e:
                    self.logger.warning('Could not convert the following InChI: '+str(chem['inchi']))
                #create a new species
                #here we want to gather the info from rpReader's rp_strc and cid_strc
                try:
                    chem_name = self.cid_strc[meta]['name']
                except KeyError:
                    chem_name = meta
                #compile as much info as you can
                #xref
                try:
                    #TODO: add the xref from the document
                    spe_xref = self.cid_xref[meta]
                except KeyError:
                    #spe_xref = {}
                    spe_xref = chem['dbref']
                #inchi
                try:
                    spe_inchi = self.cid_strc[meta]['inchi']
                except KeyError:
                    spe_inchi = chem['inchi']
                #inchikey
                try:
                    spe_inchikey = self.cid_strc[meta]['inchikey']
                except KeyError:
                    spe_inchikey =  resConv['inchikey']
                #smiles
                try:
                    spe_smiles = self.cid_strc[meta]['smiles']
                except KeyError:
                    spe_smiles = resConv['smiles']
                #pass the information to create the species
                rpsbml.createSpecies(meta,
                                     compartment_id,
                                     chem_name,
                                     spe_xref,
                                     spe_inchi,
                                     spe_inchikey,
                                     spe_smiles,
                                     species_group_id)
            #4) add the complete reactions and their annotations
            #create a new group for the measured pathway
            #need to convert the validation to step for reactions
            for stepNum in data[path_id]['steps']:
                toSend = {'left': {}, 'right': {}, 'rule_id': None, 'rule_ori_reac': None, 'rule_score': None, 'path_id': path_id, 'step': stepNum, 'sub_step': None}
                for chem in data[path_id]['steps'][stepNum]['substrates']:
                    if 'mnx' in chem['dbref']:
                        meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                        #try CHEBI
                    else:
                        self.logger.warning('Not all the species to have a MNX ID')
                        #break
                        try:
                            meta = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                            meta = 'CHEBI_'+str(meta)
                        except KeyError:
                            #TODO: need to find a better way
                            self.logger.warning('Cannot determine MNX or CHEBI entry, using random')
                            tmpDB_name = list(chem['dbref'].keys())[0]
                            meta = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                            meta = str(tmpDB_name)+'_'+str(meta) 
                    toSend['left'][meta] = 1
                for chem in data[path_id]['steps'][stepNum]['products']:
                    if 'mnx' in chem['dbref']:
                        meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                        #try CHEBI
                    else:
                        self.logger.warning('Need all the species to have a MNX ID')
                        try:
                            meta = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                            meta = 'CHEBI_'+str(meta)
                        except KeyError:
                            #TODO: need to find a better way
                            self.logger.warning('Cannot determine MNX or CHEBI entry, using random')
                            tmpDB_name = list(chem['dbref'].keys())[0]
                            meta = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                            meta = str(tmpDB_name)+'_'+str(meta) 
                    toSend['right'][meta] = 1
                        #break
                #if all are full add it
                reac_xref = {}
                if 'ec_numbers' in data[path_id]['steps'][stepNum]:
                    reac_xref['ec'] = data[path_id]['steps'][stepNum]['ec_numbers']
                if 'uniprot' in data[path_id]['steps'][stepNum]:
                    reac_xref['uniprot'] = data[path_id]['steps'][stepNum]['uniprot']
                self.logger.debug('#########################################')
                self.logger.debug(toSend)
                self.logger.debug('#########################################')
                rpsbml.createReaction(header_name+'_Step'+str(stepNum),
                                      upper_flux_bound,
                                      lower_flux_bound,
                                      toSend,
                                      compartment_id,
                                      None,
                                      reac_xref,
                                      pathway_id)
                if stepNum==1:
                    #adding the consumption of the target
                    targetStep = {'rule_id': None,
                                  'left': {},
                                  'right': {},
                                  'step': None,
                                  'sub_step': None,
                                  'path_id': None,
                                  'transformation_id': None,
                                  'rule_score': None,
                                  'rule_ori_reac': None}
                    for chem in data[path_id]['steps'][stepNum]['products']:
                        try:
                            #smallest MNX
                            meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                        except KeyError:
                            #try CHEBI
                            try:
                                meta = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                                meta = 'CHEBI_'+str(meta)
                            except KeyError:
                                self.logger.warning('Cannot determine MNX or CHEBI entry, using random')
                                tmpDB_name = list(chem['dbref'].keys())[0]
                                meta = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                                meta = str(tmpDB_name)+'_'+str(meta)
                        targetStep['left'][meta] = 1
                    rpsbml.createReaction(header_name+'_Step1_sink',
                                          upper_flux_bound,
                                          lower_flux_bound,
                                          targetStep,
                                          compartment_id)
                    rpsbml.createFluxObj('rpFBA_obj', header_name+'_Step1_sink', 1, True)
            if tmpOutputFolder:
                rpsbml.writeSBML(tmpOutputFolder)
            else:
                sbml_paths[header_name+'_Path'+str(path_id)] = rpsbml
        if tmpOutputFolder:
            return {}
        else:
            return sbml_paths


    #######################################################################
    ############################# JSON input ##############################
    #######################################################################


    #WARNING: we are not setting any restrictions on the number of steps allowed here and instead we
    #take the rule with the highest dimension. Also assume that there is a single rule at a maximal
    #dimension
    ## Function to generate an SBLM model from a JSON file
    #
    #  Read the json files of a folder describing pathways and generate an SBML file for each
    #  TODO: remove the default MNXC3 compartment ID
    #  TODO: change the ID of all species to take a normal string and not sepcial caracters
    #  WARNING: We are only using a single rule (technically with the highest diameter)
    #
    #  @param self Object pointer
    # @param colJson Dictionnary of
    #  @return rpsbml.document the SBML document
    #TODO: update this to include _hdd parsing
    #TODO: add the new group of sink species
    def jsonToSBML(self,
                   collJson,
                   upper_flux_bound=999999,
                   lower_flux_bound=0,
                   pathway_id='rp_pathway',
                   compartment_id='MNXC3',
                   species_group_id='central_species',
                   pubchem_search=False):
        """Function to generate an SBLM model from a JSON file (RetroPathRL)

        :param collJson: The input JSON file
        :param upper_flux_bound: The default upper flux bound (Default: 999999)
        :param lower_flux_bound: The default lower flux bound (Default: 0)
        :param pathway_id: The Groups heterologous pathway id (Default: rp_pathway)
        :param compartment_id: The compartment SBML id (Default: MNXC3)
        :param species_group_id: The Groups id of the central species (Default: central_species)
        :param pubchem_search: Use the pubchem database to search for missing cross reference (Default: False)

        :type collJson: str
        :type upper_flux_bound: int
        :type lower_flux_bound: int
        :type pathway_id: str
        :type compartment_id: str
        :type species_group_id: str
        :type pubchem_search: str

        :rtype: dict
        :return: The dictionary of the pathways
        """
        #global parameters used for all parameters
        pathNum = 1
        rp_paths = {}
        reac_smiles = {}
        reac_ecs = {}
        species_list = {}
        reactions_list = {}
        source_cid = {}
        source_stochio = {}
        cid_inchikey = {}
        sink_species = {}
        ############################################
        ############## gather the data #############
        ############################################
        for json_dict in collJson:
            ########### construct rp_paths ########
            reac_smiles[pathNum] = {}
            reac_ecs[pathNum] = {}
            species_list[pathNum] = {}
            reactions_list[pathNum] = {}
            cid_inchikey[pathNum] = {}
            sink_species[pathNum] = {}
            stochio = {}
            inchikey_cid = {}
            source_species = []
            skip_pathway = False
            for node in collJson[json_dict]['elements']['nodes']:
                ##### compounds ####
                if node['data']['type']=='compound':
                    cid_inchikey[pathNum][node['data']['id'].replace('-', '')] = node['data']['id']
                    species_list[pathNum][node['data']['id'].replace('-', '')] = {'inchi': node['data']['InChI'],
                                    'inchikey': node['data']['id'],
                                    'smiles': node['data']['SMILES']}
                    if int(node['data']['isSource'])==1:
                        #TODO: there should always be only one source, to check
                        source_species.append(node['data']['id'].replace('-', ''))
                        source_cid[pathNum] = node['data']['id'].replace('-', '')
                    if int(node['data']['inSink'])==1:
                        sink_species[pathNum][node['data']['id'].replace('-', '')] = node['data']['id']
                ###### reactions ######
                elif node['data']['type']=='reaction':
                    #NOTE: pick the rule with the highest diameter
                    r_id = sorted(node['data']['Rule ID'], key=lambda x: int(x.split('-')[-2]), reverse=True)[0]
                    reactions_list[pathNum][node['data']['id']] = {'rule_id': r_id,
                                                                   'rule_ori_reac': None,
                                                                   'right': {},
                                                                   'left': {},
                                                                   'path_id': pathNum,
                                                                   'step': None,
                                                                   'sub_step': None,
                                                                   'transformation_id': None,
                                                                   'rule_score': node['data']['Score']}
                    reac_smiles[pathNum][r_id] = node['data']['Reaction SMILES']
                    reac_ecs[pathNum][r_id] = list(filter(None, [i for i in node['data']['EC number']]))
                    stochio[node['data']['id']] = {}
                    for i in node['data']['Stoechiometry']:
                        stochio[node['data']['id']][i.replace('-', '')] = node['data']['Stoechiometry'][i]
            ############ make main pathway ###########
            main_reac = {}
            for reaction_node in collJson[json_dict]['elements']['edges']:
                if not len(reaction_node['data']['source'].split('-'))==3:
                    if not reaction_node['data']['source'] in reactions_list[pathNum]:
                        self.logger.error('The following reaction was not found in the JSON elements: '+str(reaction_node['data']['source']))
                        skip_pathway = True
                        break
                    else:
                        rid = reaction_node['data']['source']
                        try:
                            cid = inchikey_cid[reaction_node['data']['target'].replace('-', '')]
                        except KeyError:
                            cid = reaction_node['data']['target'].replace('-', '')
                        try:
                            reactions_list[pathNum][rid]['left'][cid] = stochio[rid][cid]
                        except KeyError:
                            reactions_list[pathNum][rid]['left'][cid] = 1.0
                            self.logger.warning('The cid ('+str(cid)+') has not been detected by stochio. Setting to 1.0')
                if not len(reaction_node['data']['target'].split('-'))==3:
                    if not reaction_node['data']['target'] in reactions_list[pathNum]:
                        self.logger.error('The following reaction was not found in the JSON elements: '+str(reaction_node['data']['source']))
                        skip_pathway = True
                        break
                    else:
                        rid = reaction_node['data']['target']
                        try:
                            cid = inchikey_cid[reaction_node['data']['source'].replace('-', '')]
                        except KeyError:
                            cid = reaction_node['data']['source'].replace('-', '')
                        try:
                            reactions_list[pathNum][rid]['right'][cid] = stochio[rid][cid]
                        except KeyError:
                            reactions_list[pathNum][rid]['right'][cid] = 1.0
                            self.logger.warning('The cid ('+str(cid)+') has not been detected by stochio. Setting to 1.0')
            ################# calculate the steps associated with the reactions_list ######
            #find the source in the LAST reaction in the pathway
            #NOTE: this assumes that the source is contained in a single final reaction and nowhere else
            last_rid = None
            step_num = 1
            found_rid = []
            toFind_rid = list(reactions_list[pathNum].keys())
            for rid in reactions_list[pathNum]:
                if all([True if i in source_species else False for i in reactions_list[pathNum][rid]['right']]):
                    reactions_list[pathNum][rid]['step'] = step_num
                    #step_num -= 1
                    step_num += 1
                    last_rid = rid
                    try:
                        source_stochio[pathNum] = stochio[rid][source_cid[pathNum]]
                    except KeyError:
                        source_stochio[pathNum] = 1.0
                    found_rid.append(rid)
                    toFind_rid.remove(rid)
                    break
            for rid in toFind_rid[:]:
                if all([True if i in reactions_list[pathNum][last_rid]['left'] else False for i in reactions_list[pathNum][rid]['right']]):
                    reactions_list[pathNum][rid]['step'] = step_num
                    #step_num -= 1
                    step_num += 1
                    last_rid = rid
                    found_rid.append(rid)
                    toFind_rid.remove(rid)
            if not toFind_rid==[]:
                self.logger.error('There are reactions unaccounted for: '+str(toFind_rid))
                skip_pathway = True
                break
            ############# find all the alternative reactions associated with a reaction rule ###
            if not skip_pathway:
                rp_paths[pathNum] = {}
                for rid in reactions_list[pathNum]:
                    rp_paths[pathNum][reactions_list[pathNum][rid]['step']] = {}
                    sub_step = 1
                    for reac_id in self.rr_reactions[reactions_list[pathNum][rid]['rule_id']]:
                        tmpReac = copy.deepcopy(reactions_list[pathNum][rid])
                        tmpReac['mnxr'] = reac_id
                        tmpReac['sub_step'] = sub_step
                        rp_paths[pathNum][reactions_list[pathNum][rid]['step']][sub_step] = tmpReac
                        sub_step += 1
            else:
                self.logger.warning('Skipping pathway '+str(pathNum))
            pathNum += 1
        ######################################
        ########### create the SBML's ########
        ######################################
        try:
            mnxc = self.xref_comp[compartment_id]
        except KeyError:
            self.logger.error('Could not Xref compartment_id ('+str(compartment_id)+')')
            return False
        sbml_paths = {}
        for pathNum in rp_paths:
            #first level is the list of lists of sub_steps
            #second is itertools all possible combinations using product
            altPathNum = 1
            for comb_path in list(itertools.product(*[[(i,y) for y in rp_paths[pathNum][i]] for i in rp_paths[pathNum]])):
                steps = []
                for i, y in comb_path:
                    steps.append(rp_paths[pathNum][i][y])
                path_id = steps[0]['path_id']
                rpsbml = rpSBML.rpSBML('rp_'+str(path_id)+'_'+str(altPathNum))
                #1) create a generic Model, ie the structure and unit definitions that we will use the most
                ##### TODO: give the user more control over a generic model creation:
                #   -> special attention to the compartment
                rpsbml.genericModel('RetroPath_Pathway_'+str(path_id)+'_'+str(altPathNum),
                                    'RP_model_'+str(path_id)+'_'+str(altPathNum),
                                    self.comp_xref[mnxc],
                                    compartment_id,
                                    upper_flux_bound,
                                    lower_flux_bound)
                #2) create the pathway (groups)
                rpsbml.createPathway(pathway_id)
                rpsbml.createPathway(species_group_id)
                #3) find all the unique species and add them to the model
                ###################################################
                ############## SPECIES ############################
                ###################################################
                meta_to_cid = {}
                for meta in species_list[pathNum]:
                    #### beofre adding it to the model check to see if you can recover some MNXM from inchikey
                    #NOTE: only for the sink species do we try to convert to MNXM
                    if meta in list(sink_species[pathNum].keys()):
                        try:
                            #take the smallest MNX, usually the best TODO: review this
                            cid = sorted(self.inchikey_cid[sink_species[pathNum][meta]], key=lambda x: int(x[4:]))[0]
                            meta_to_cid[meta] = cid
                        except KeyError:
                            self.logger.warning('Cannot find sink compound: '+str(meta))
                            continue
                            #return False
                    else:
                        cid = meta
                    # retreive the name of the molecule
                    #here we want to gather the info from rpReader's rp_strc and cid_strc
                    try:
                        chem_name = self.cid_strc[meta]['name']
                    except KeyError:
                        chem_name = None
                    #compile as much info as you can
                    #xref
                    try:
                        spe_xref = self.cid_xref[meta]
                    except KeyError:
                        spe_xref = {}
                    ###### Try to recover the structures ####
                    spe_smiles = None
                    spe_inchi = None
                    spe_inchikey = None
                    pubchem_smiles = None
                    pubchem_inchi = None
                    pubchem_inchikey = None
                    pubchem_xref = {}
                    #inchi
                    try:
                        spe_inchi = rp_strc[meta]['inchi']
                        if not spe_xref and pubchem_search:
                            try:
                                if not self.pubchem_inchi[spe_inchi]=={}:
                                    pubchem_inchi = self.pubchem_inchi[spe_inchi]['inchi']
                                    pubchem_inchikey = self.pubchem_inchi[spe_inchi]['inchikey']
                                    pubchem_smiles = self.pubchem_inchi[spe_inchi]['smiles']
                                    pubchem_xref = self.pubchem_inchi[spe_inchi]['xref'] 
                            except KeyError:
                                pubres = self._pubchemStrctSearch(spe_inchi, 'inchi')
                                if not chem_name:
                                    chem_name = pubres['name']
                                if 'chebi' in pubres['xref']:
                                    try:
                                        #WARNING: taking the first one. Better to take the smallest?
                                        spe_xref = self.cid_xref[self.chebi_cid[pubres['xref']['chebi'][0]]]
                                    except KeyError:
                                        pass
                                if not pubchem_xref:
                                    pubchem_xref = pubres['xref']
                                if not pubchem_inchikey:
                                    pubchem_inchikey = pubres['inchikey']
                                if not pubchem_smiles:
                                    pubchem_smiles = pubres['smiles']
                    except KeyError:
                        self.logger.warning('Bad results from pubchem results')
                        pass
                    #inchikey
                    try:
                        spe_inchikey = rp_strc[meta]['inchikey']
                        if not spe_xref and pubchem_search:
                            pubres = self._pubchemStrctSearch(spe_inchikey, 'inchikey')
                            if not chem_name:
                                chem_name = pubres['name']
                            if 'chebi' in pubres['xref']:
                                try:
                                    spe_xref = self.cid_xref[self.chebi_cid[pubres['xref']['chebi'][0]]]
                                except KeyError:
                                    pass
                            if not pubchem_xref:
                                pubchem_xref = pubres['xref']
                            if not pubchem_inchi:
                                pubchem_inchi = pubres['inchi']
                            if not pubchem_smiles:
                                pubchem_smiles = pubres['smiles']
                    except KeyError:
                        self.logger.warning('Bad results from pubchem results')
                        pass
                    #smiles
                    try:
                        spe_smiles = rp_strc[meta]['smiles']
                        if not spe_xref and pubchem_search:
                            pubres = self._pubchemStrctSearch(spe_smiles, 'smiles')
                            if not chem_name:
                                chem_name = pubres['name']
                            if 'chebi' in pubres['xref']:
                                try:
                                    spe_xref = self.cid_xref[self.chebi_cid[pubres['xref']['chebi'][0]]]
                                except KeyError:
                                    pass        
                            if not pubchem_xref:
                                pubchem_xref = pubres['xref']
                            if not pubchem_inchi:
                                pubchem_inchi = pubres['inchi']
                            if not pubchem_inchikey:
                                pubchem_inchikey = pubres['inchikey']
                    except KeyError:
                        self.logger.warning('Bad results from pubchem results')
                        pass
                    if not spe_inchi:
                        spe_inchi = pubchem_inchi
                    if not spe_inchikey:
                        spe_inchikey = pubchem_inchikey
                    if not spe_smiles:
                        spe_smiles = pubchem_smiles
                    if not spe_xref:
                        spe_xref = pubchem_xref
                    #pass the information to create the species
                    rpsbml.createSpecies(meta,
                            compartment_id,
                            chem_name,
                            spe_xref,
                            spe_inchi,
                            spe_inchikey,
                            spe_smiles,
                            species_group_id)
                #4) add the reactions and their annotations
                for step in steps:
                    #add the substep to the model
                    step['sub_step'] = altPathNum
                    #### try to replace the sink compounds inchikeys with mnxm
                    for direc in ['right', 'left']:
                        step_mnxm = {}
                        for meta in step[direc]:
                            try:
                                step_mnxm[meta_to_cid[meta]] = step[direc][meta]
                            except KeyError:
                                step_mnxm[meta] = step[direc][meta]
                        step[direc] = step_mnxm
                    rpsbml.createReaction('RP'+str(step['step']),
                            upper_flux_bound,
                            lower_flux_bound,
                            step,
                            compartment_id,
                            reac_smiles[pathNum][step['rule_id']],
                            {'ec': reac_ecs[pathNum][step['rule_id']]},
                            pathway_id)
                #5) adding the consumption of the target
                targetStep = {'rule_id': None,
                        'left': {source_cid[pathNum]: source_stochio[pathNum]},
                        'right': {},
                        'step': None,
                        'sub_step': None,
                        'path_id': None,
                        'transformation_id': None,
                        'rule_score': None,
                        'rule_ori_reac': None}
                rpsbml.createReaction('RP1_sink',
                        upper_flux_bound,
                        lower_flux_bound,
                        targetStep,
                        compartment_id)
                #6) Optional?? Add the flux objectives. Could be in another place, TBD
                rpsbml.createFluxObj('rpFBA_obj', 'RP1_sink', 1, True)
                sbml_paths['rp_'+str(step['path_id'])+'_'+str(altPathNum)] = rpsbml
                altPathNum += 1



    '''
    ######################################################
    ################## string to sbml ####################
    ######################################################
    ##
    # react_string: '1 MNX:MNXM181 + 1 MNX:MNXM4 => 2 MNX:MNXM1 + 1 MNX:MNXM11441'
    # ec: []
    def reacStr_to_sbml(self,
                        reac_sctring,
                        ec=[])
        #[{'reactants': [{'inchi': 'ajsjsjsjs', 'db_name': 'mnx', 'name': 'species1', 'id': 'MNXM181'}], 'products': [{'inchi': 'jdjdjdjdj', 'db_name': 'mnx', 'name': 'product1', 'id': 'MNXM1'}], 'ec': [{'id': '1.1.1.1'}]}]
        reac_species = [{'stoichio': i.split(' ')[0],
                         'inchi': '',
                         'name': i.split(' ')[0].split(':')[1],
                         'db_name': i.split(' ')[0].split(':')[0].lower(),
                         'id': i.split(' ')[0].split(':')[1]}  for i in reac_sctring.split('=>')[0].split('+')]
        reac_products = [{'stoichio': i.split(' ')[1],
                          'inchi': '',
                          'name': i.split(' ')[1].split(':')[1],
                          'db_name': i.split(' ')[1].split(':')[0].lower(),
                          'id': i.split(' ')[1].split(':')[1]}  for i in reac_sctring.split('=>')[0].split('+')]
        reac = [{'reactants': reac_species,
                 'products': reac_products,
                 'ec': [{'id': i} for i in ec]}]
	########### create CSV from input ##########
	with tempfile.TemporaryDirectory() as tmpOutputFolder:
		with open(tmpOutputFolder+'/tmp_input.tsv', 'w') as infi:
			csvfi = csv.writer(infi, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
			header = ['pathway_ID',
                                  'target_name',
                                  'target_structure',
                                  'step',
                                  'substrate_name',
                                  'substrate_dbref',
                                  'substrate_structure',
                                  'product_name',
                                  'product_dbref',
                                  'product_structure',
                                  'EC_number',
                                  'enzyme_identifier',
                                  'enzyme_name',
                                  'organism',
                                  'yield',
                                  'comments',
                                  'pictures',
                                  'pdf',
                                  'reference',
                                  'estimated time',
                                  'growth media']
			csvfi.writerow(header)
			first_line = ['1',
                                      'void',
                                      'test',
                                      '0',
                                      '',
                                      '',
                                      '',
                                      '',
                                      '',
                                      '',
                                      '',
                                      '',
                                      '',
                                      '',
                                      '',
                                      '',
                                      '',
                                      '',
                                      '',
                                      '',
                                      '']
			csvfi.writerow(first_line)
			reac_count = len(reac)
			for reaction in reac:
				reac = ';'.join([i['id'] for i in reaction['reactants']])
				to_write = ['1',
					input_dict['target_name'],
					input_dict['target_inchi'],
					str(reac_count),
					';'.join([i['name'] for i in reaction['reactants']]),
					';'.join([i['db_name']+':'+i['id'] for i in reaction['reactants']]),
					';'.join([i['inchi'] for i in reaction['reactants']]),
					';'.join([i['name'] for i in reaction['products']]),
					';'.join([i['db_name']+':'+i['id'] for i in reaction['products']]),
					';'.join([i['inchi'] for i in reaction['products']]),
					';'.join([i['id'] for i in reaction['ec']]),
					'',
					'',
					'',
					'',
					'',
					'',
					'',
					'',
					'',
					'']
				csvfi.writerow(to_write)
				reac_count -= 1
		############## create SBML from CSV #####
		### rpReader #####
		rpreader = rpReader.rpReader()
		rpcache = rpToolCache.rpToolCache()
		rpreader.deprecatedCID_cid = rpcache.deprecatedCID_cid
		rpreader.deprecatedRID_rid = rpcache.deprecatedRID_rid
		rpreader.cid_strc = rpcache.cid_strc
		rpreader.inchikey_cid = rpcache.inchikey_cid
		rpreader.rr_reactions = rpcache.rr_reactions
		rpreader.cid_xref = rpcache.cid_xref
		rpreader.comp_xref = rpcache.comp_xref
		rpreader.xref_compID = rpcache.xref_compID
		##################
		#measured_pathway = parseValidation(tmpOutputFolder+'/tmp_input.csv')
		measured_sbml_paths = rpreader.validationToSBML(tmpOutputFolder+'/tmp_input.tsv',
                                                                tmpOutputFolder+'/')
		#should be only one
		measured_sbml = glob.glob(tmpOutputFolder+'/*.sbml')[0]
    '''
