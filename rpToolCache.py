from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs
import csv
import logging
import os
import pickle
import gzip
import urllib.request
from rpCache import rpCache

## Error function for the convertion of structures
#
class Error(Exception):
    pass


## Error function for the convertion of structures
#
class DepictionError(Error):
    def __init__(self, message):
        #self.expression = expression
        self.message = message


#######################################################
################### rpCache  ##########################
#######################################################

## Class to generate the cache
#
# Contains all the functions that parse different files, used to calculate the thermodynamics and the FBA of the
#the other steps. These should be called only when the files have changes
class rpToolCache(rpCache):
    ## Cache constructor
    #
    # @param self The object pointer
    # @param inputPath The path to the folder that contains all the input/output files required
    def __init__(self):

        super().__init__()

        #given by Thomas
        self.inchikey_mnxm = None
        self.compXref = None
        self.nameCompXref = None
        if not self._loadCache():
            raise KeyError



    ## Private function to fetch the required data, parse them and generate the pickle
    #
    #  Opens the previously generated cache to the object memory
    #
    # @param The oject pointer
    # @return Boolean detemining the success of the function or not
    def _loadCache(self, fetchInputFiles=False):

        dirname = os.path.dirname(os.path.abspath( __file__ ))
        # comp_xref
        if not os.path.isfile(dirname+'/input_cache/comp_xref.tsv') or fetchInputFiles:
            urllib.request.urlretrieve('https://www.metanetx.org/cgi-bin/mnxget/mnxref/comp_xref.tsv',
                                       dirname+'/input_cache/comp_xref.tsv')

        ###################### Populate the cache #################################
        if not os.path.isfile(dirname+'/cache/inchikey_mnxm.pickle.gz'):
            pickle.dump(self.mnx_strc(dirname+'/input_cache/rr_compounds.tsv',
                                      dirname+'/input_cache/chem_prop.tsv'),
                        gzip.open(dirname+'/cache/inchikey_mnxm.pickle.gz','wb'))
        self.inchikey_mnxm = pickle.load(gzip.open(dirname+'/cache/inchikey_mnxm.pickle.gz', 'rb'))

        if not os.path.isfile(dirname+'/cache/inchikey_mnxm.pickle.gz'):
            inchikey_mnxm = {}
            for mnxm in self.mnxm_strc:
                if not self.mnxm_strc[mnxm]['inchikey'] in inchikey_mnxm:
                    inchikey_mnxm[self.mnxm_strc[mnxm]['inchikey']] = []
                inchikey_mnxm[self.mnxm_strc[mnxm]['inchikey']].append(mnxm)
            pickle.dump(inchikey_mnxm, gzip.open(dirname+'/cache/inchikey_mnxm.pickle.gz','wb'))

        if not os.path.isfile(dirname+'/cache/compXref.pickle.gz') or not os.path.isfile(dirname+'/cache/nameCompXref.pickle.gz'):
            name_pubDB_xref, compName_mnxc = self.mnx_compXref(dirname+'/input_cache/comp_xref.tsv')
            pickle.dump(name_pubDB_xref, gzip.open(dirname+'/cache/compXref.pickle.gz','wb'))
            pickle.dump(compName_mnxc, gzip.open(dirname+'/cache/nameCompXref.pickle.gz','wb'))
        self.compXref = pickle.load(gzip.open(dirname+'/cache/compXref.pickle.gz', 'rb'))

        self.nameCompXref = pickle.load(gzip.open(dirname+'/cache/nameCompXref.pickle.gz', 'rb'))

        return super()._loadCache(fetchInputFiles)


    #######################################################
    ################### PRIVATE FUNCTION ##################
    #######################################################



    ## Function to parse the compXref.tsv file of MetanetX
    #
    #  Generate a dictionnary of compartments id's (MNX) to other database id's
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def mnx_compXref(self, compXref_path):
        name_pubDB_xref = {}
        compName_mnxc = {}
        try:
            with open(compXref_path) as f:
                c = csv.reader(f, delimiter='\t')
                #not_recognised = []
                for row in c:
                    #cid = row[0].split(':')
                    if not row[0][0]=='#':
                        #collect the info
                        mnxc = row[1]
                        if len(row[0].split(':'))==1:
                            dbName = 'mnx'
                            dbCompId = row[0]
                        else:
                            dbName = row[0].split(':')[0]
                            dbCompId = ''.join(row[0].split(':')[1:])
                            dbCompId = dbCompId.lower()
                        if dbName=='deprecated':
                            dbName = 'mnx'
                        #create the dicts
                        if not mnxc in name_pubDB_xref:
                            name_pubDB_xref[mnxc] = {}
                        if not dbName in name_pubDB_xref[mnxc]:
                            name_pubDB_xref[mnxc][dbName] = []
                        if not dbCompId in name_pubDB_xref[mnxc][dbName]:
                            name_pubDB_xref[mnxc][dbName].append(dbCompId)
                        #create the reverse dict
                        if not dbCompId in compName_mnxc:
                            compName_mnxc[dbCompId] = mnxc
        except FileNotFoundError:
            self.logger.error('compXref file not found')
            return {}
        return name_pubDB_xref, compName_mnxc


if __name__ == "__main__":
    rpcache = rpToolCache()
