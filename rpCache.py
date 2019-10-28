from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs
import csv

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
class rpCache:
    ## Cache constructor
    #
    # @param self The object pointer
    # @param inputPath The path to the folder that contains all the input/output files required
    def __init__(self):
        #given by Thomas
        self.logger = logging.getLogger(__name__)
        self.logger.info('Started instance of rpCache')
        self.convertMNXM = {'MNXM162231': 'MNXM6',
                'MNXM84': 'MNXM15',
                'MNXM96410': 'MNXM14',
                'MNXM114062': 'MNXM3',
                'MNXM145523': 'MNXM57',
                'MNXM57425': 'MNXM9',
                'MNXM137': 'MNXM588022'}
        self.deprecatedMNXM_mnxm = {}
        self.deprecatedMNXR_mnxr = {}
        self.mnxm_strc = None
        self.inchikey_mnxm = None
        self.rr_reactions = None
        self.chemXref = None
        self.compXref = None
        self.nameCompXref = None



    ## Private function to fetch the required data, parse them and generate the pickle
    #
    #  Opens the previously generated cache to the object memory
    #
    # @param The oject pointer
    # @return Boolean detemining the success of the function or not
    def _loadCache(self, fetchInputFiles=False):
        dirname = os.path.dirname(os.path.abspath( __file__ ))
        #################### make the local folders ############################
        # input_cache
        if not os.path.isdir(dirname+'/input_cache'):
            os.mkdir(dirname+'/input_cache')
        # cache
        if not os.path.isdir(dirname+'/cache'):
            os.mkdir(dirname+'/cache')
        ###################### Fetch the files if necessary ######################
        #chem_xref
        if not os.path.isfile(dirname+'/input_cache/chem_xref.tsv') or fetchInputFiles:
            urllib.request.urlretrieve('https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv', 
                    dirname+'/input_cache/chem_xref.tsv')
        #reac_xref
        if not os.path.isfile(dirname+'/input_cache/reac_xref.tsv') or fetchInputFiles:
            urllib.request.urlretrieve('https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv', 
                    dirname+'/input_cache/reac_xref.tsv')
        # rules_rall.tsv
        if not os.path.isfile(dirname+'/input_cache/rules_rall.tsv') or fetchInputFiles:
            urllib.request.urlretrieve('TODO', 
                    dirname+'/input_cache/rules_rall.tsv')
            '''
            tf = tarfile.open(dirname+'/input_cache/retrorules_preparsed.tar.gz')
            tf.extractall(path=dirname+'/input_cache/')
            tf.close()
            os.rename(dirname+'/retrorules_preparsed/retrorules-rr01_flat_all.csv', 
                    dirname+'/input_cache/rules_rall.tsv')
            '''
        # rr_compounds.tsv
        #TODO: need to add this file to the git or another location
        if not os.path.isfile(dirname+'/input_cache/rr_compounds.tsv') or fetchInputFiles:
            urllib.request.urlretrieve(
                    'TODO', 
                    dirname+'/input_cache/rr_compounds.tsv')
            '''
            tf = tarfile.open(dirname+'/input_cache/retrorules_preparsed.tar.xz')
            tf.extractall(path=dirname+'/input_cache/')
            tf.close()
            '''
        # chem_prop.tsv
        if not os.path.isfile(dirname+'/input_cache/chem_prop.tsv') or fetchInputFiles:
            urllib.request.urlretrieve('https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv', 
                    dirname+'/input_cache/chem_prop.tsv')
        # comp_xref
        if not os.path.isfile(dirname+'/input_cache/comp_xref.tsv') or fetchInputFiles:
            urllib.request.urlretrieve('https://www.metanetx.org/cgi-bin/mnxget/mnxref/comp_xref.tsv', 
                    dirname+'/input_cache/comp_xref.tsv')
        ###################### Populate the cache #################################
        rpcache = rpCache()
        if not os.path.isfile(dirname+'/cache/deprecatedMNXM_mnxm.pickle'):
            rpcache.deprecatedMNXM(dirname+'/input_cache/chem_xref.tsv')
            pickle.dump(rpcache.deprecatedMNXM_mnxm, open(dirname+'/cache/deprecatedMNXM_mnxm.pickle', 'wb'))
        self.deprecatedMNXM_mnxm = pickle.load(open(dirname+'/cache/deprecatedMNXM_mnxm.pickle', 'rb'))
        if not os.path.isfile(dirname+'/cache/deprecatedMNXR_mnxr.pickle'):
            rpcache.deprecatedMNXR(dirname+'/input_cache/reac_xref.tsv')
            pickle.dump(rpcache.deprecatedMNXR_mnxr, open(dirname+'/cache/deprecatedMNXR_mnxr.pickle', 'wb'))
        self.deprecatedMNXR_mnxr = pickle.load(open(dirname+'/cache/deprecatedMNXR_mnxr.pickle', 'rb'))
        if not os.path.isfile(dirname+'/cache/mnxm_strc.pickle.gz') or not os.path.isfile(dirname+'/cache/inchikey_mnxm.pickle.gz'):
            mnxm_strc, inchikey_mnxm = rpcache.mnx_strc(dirname+'/input_cache/rr_compounds.tsv',
                                                        dirname+'/input_cache/chem_prop.tsv')
            pickle.dump(mnxm_strc, gzip.open(dirname+'/cache/mnxm_strc.pickle.gz','wb'))
            pickle.dump(inchikey_mnxm, gzip.open(dirname+'/cache/inchikey_mnxm.pickle.gz','wb'))
        self.inchikey_mnxm = pickle.load(gzip.open(dirname+'/cache/inchikey_mnxm.pickle.gz', 'rb'))
        self.mnxm_strc = pickle.load(gzip.open(dirname+'/cache/mnxm_strc.pickle.gz', 'rb'))
        if not os.path.isfile(dirname+'/cache/inchikey_mnxm.pickle.gz'):
            inchikey_mnxm = {}
            for mnxm in self.mnxm_strc:
                if not self.mnxm_strc[mnxm]['inchikey'] in inchikey_mnxm:
                    inchikey_mnxm[self.mnxm_strc[mnxm]['inchikey']] = []
                inchikey_mnxm[self.mnxm_strc[mnxm]['inchikey']].append(mnxm)
            pickle.dump(inchikey_mnxm, gzip.open(dirname+'/cache/inchikey_mnxm.pickle.gz','wb'))
        if not os.path.isfile(dirname+'/cache/rr_reactions.pickle'):
            pickle.dump(
                    rpcache.retro_reactions(dirname+'/input_cache/rules_rall.tsv'), 
                    open(dirname+'/cache/rr_reactions.pickle', 'wb'))
        self.rr_reactions = pickle.load(open(dirname+'/cache/rr_reactions.pickle', 'rb'))
        if not os.path.isfile(dirname+'/cache/chemXref.pickle.gz'):
            pickle.dump(rpcache.mnx_chemXref(dirname+'/input_cache/chem_xref.tsv'), 
                    gzip.open(dirname+'/cache/chemXref.pickle.gz','wb'))
        self.chemXref = pickle.load(gzip.open(dirname+'/cache/chemXref.pickle.gz', 'rb'))
        if not os.path.isfile(dirname+'/cache/compXref.pickle.gz') or not os.path.isfile(dirname+'/cache/nameCompXref.pickle.gz'):
            name_pubDB_xref, compName_mnxc = rpcache.mnx_compXref(dirname+'/input_cache/comp_xref.tsv')
            pickle.dump(name_pubDB_xref, gzip.open(dirname+'/cache/compXref.pickle.gz','wb'))
            pickle.dump(compName_mnxc, gzip.open(dirname+'/cache/nameCompXref.pickle.gz','wb'))
        self.compXref = pickle.load(gzip.open(dirname+'/cache/compXref.pickle.gz', 'rb'))
        self.nameCompXref = pickle.load(gzip.open(dirname+'/cache/nameCompXref.pickle.gz', 'rb'))
        rpcache = None
        return True


    #######################################################
    ################### PRIVATE FUNCTION ##################
    #######################################################


    ## Convert chemical depiction to others type of depictions
    #
    # Usage example:
    # - convert_depiction(idepic='CCO', otype={'inchi', 'smiles', 'inchikey'})
    # - convert_depiction(idepic='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', itype='inchi', otype={'inchi', 'smiles', 'inchikey'})
    #  @param self The object pointer
    #  @param idepic String depiction to be converted, str
    #  @param itype type of depiction provided as input, str
    #  @param otype types of depiction to be generated, {"", "", ..}
    #  @return odepic generated depictions, {"otype1": "odepic1", ..}
    def _convert_depiction(self, idepic, itype='smiles', otype={'inchikey'}):
        # Import (if needed)
        if itype == 'smiles':
            rdmol = MolFromSmiles(idepic, sanitize=True)
        elif itype == 'inchi':
            rdmol = MolFromInchi(idepic, sanitize=True)
        else:
            raise NotImplementedError('"{}" is not a valid input type'.format(itype))
        if rdmol is None:  # Check imprt
            raise DepictionError('Import error from depiction "{}" of type "{}"'.format(idepic, itype))
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


    #######################################################
    #######################################################
    #######################################################


    ## Function to create a dictionnary of old to new chemical id's
    #
    #  Generate a one-to-one dictionnary of old id's to new ones. Private function
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    def _checkMNXMdeprecated(self, mnxm):
        try:
            return self.deprecatedMNXM_mnxm[mnxm]
        except KeyError:
            return mnxm


    ## Function to create a dictionnary of old to new reaction id's
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    def _checkMNXRdeprecated(self, mnxr):
        try:
            return self.deprecatedMNXR_mnxr[mnxr]
        except KeyError:
            return mnxr


    #[TODO] merge the two functions
    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers to make sure that we always use the freshest id's.
    # This can include more than one old id per new one and thus returns a dictionnary. Private function
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return Dictionnary of identifiers
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def deprecatedMNXM(self, chem_xref_path):
        self.deprecatedMNXM_mnxm = {}
        with open(chem_xref_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = row[0].split(':')
                    if mnx[0]=='deprecated':
                        self.deprecatedMNXM_mnxm[mnx[1]] = row[1]
            self.deprecatedMNXM_mnxm.update(self.convertMNXM)
            self.deprecatedMNXM_mnxm['MNXM01'] = 'MNXM1'


    ## Function to parse the reac_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers to make sure that we always use the freshest id's.
    # This can include more than one old id per new one and thus returns a dictionnary. Private function
    #
    #  @param self Object pointer
    #  @param reac_xref_path Input file path
    #  @return Dictionnary of identifiers
    def deprecatedMNXR(self, reac_xref_path):
        self.deprecatedMNXMR_mnxr = {}
        with open(reac_xref_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = row[0].split(':')
                    if mnx[0]=='deprecated':
                        self.deprecatedMNXR_mnxr[mnx[1]] = row[1]


    ## Function to parse the chemp_prop.tsv file from MetanetX and compounds.tsv from RetroRules. Uses the InchIkey as key to the dictionnary
    #
    #  Generate a dictionnary gaving the formula, smiles, inchi and inchikey for the components
    #
    #  @param self Object pointer
    #  @param chem_prop_path Input file path
    #  @return mnxm_strc Dictionnary of formula, smiles, inchi and inchikey
    def mnx_strc(self, rr_compounds_path, chem_prop_path):
        mnxm_strc = {}
        for row in csv.DictReader(open(rr_compounds_path), delimiter='\t'):
            tmp = {'forumla':  None,
                    'smiles': None,
                    'inchi': row['inchi'],
                    'inchikey': None,
                    'mnxm': self._checkMNXMdeprecated(row['cid']),
                    'name': None}
            try:
                resConv = self._convert_depiction(idepic=tmp['inchi'], itype='inchi', otype={'smiles','inchikey'})
                for i in resConv:
                    tmp[i] = resConv[i]
            except DepictionError as e:
                self.logger.warning('Could not convert some of the structures: '+str(tmp))
                self.logger.warning(e)
            mnxm_strc[tmp['mnxm']] = tmp
        with open(chem_prop_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnxm = self._checkMNXMdeprecated(row[0])
                    tmp = {'forumla':  row[2],
                            'smiles': row[6],
                            'inchi': row[5],
                            'inchikey': row[8],
                            'mnxm': mnxm,
                            'name': row[1]}
                    for i in tmp:
                        if tmp[i]=='' or tmp[i]=='NA':
                            tmp[i] = None
                    if mnxm in mnxm_strc:
                        mnxm_strc[mnxm]['forumla'] = row[2]
                        mnxm_strc[mnxm]['name'] = row[1]
                        if not mnxm_strc[mnxm]['smiles'] and tmp['smiles']:
                            mnxm_strc[mnxm]['smiles'] = tmp['smiles']
                        if not mnxm_strc[mnxm]['inchikey'] and tmp['inchikey']:
                            mnxm_strc[mnxm]['inchikey'] = tmp['inchikey']
                    else:
                        #check to see if the inchikey is valid or not
                        otype = set({})
                        if not tmp['inchikey']:
                            otype.add('inchikey')
                        if not tmp['smiles']:
                            otype.add('smiles')
                        if not tmp['inchi']:
                            otype.add('inchi')
                        itype = ''
                        if tmp['inchi']:
                            itype = 'inchi'
                        elif tmp['smiles']:
                            itype = 'smiles'
                        else:
                            self.logger.warning('No valid entry for the convert_depiction function')
                            continue
                        try:
                            resConv = self._convert_depiction(idepic=tmp[itype], itype=itype, otype=otype)
                            for i in resConv:
                                tmp[i] = resConv[i]
                        except DepictionError as e:
                            self.logger.warning('Could not convert some of the structures: '+str(tmp))
                            self.logger.warning(e)
                        mnxm_strc[tmp['mnxm']] = tmp
            #inchikey_mnam
            inchikey_mnxm = {}
            for mnxm in mnxm_strc:
                if not mnxm_strc[mnxm]['inchikey'] in inchikey_mnxm:
                    inchikey_mnxm[mnxm_strc[mnxm]['inchikey']] = []
                inchikey_mnxm[mnxm_strc[mnxm]['inchikey']].append(mnxm)
        return mnxm_strc, inchikey_mnxm


    ## Function to parse the rules_rall.tsv from RetroRules
    #
    #  Extract from the reactions rules the ruleID, the reactionID, the direction of the rule directed to the origin reaction
    #
    #  @param self The object pointer.
    #  @param path The input file path.
    #  @return rule Dictionnary describing each reaction rule
    def retro_reactions(self, rules_rall_path):
        try:
            #with open(rules_rall_path, 'r') as f:
            #    reader = csv.reader(f, delimiter = '\t')
            #    next(reader)
            #    rule = {}
            #    for row in reader:
            rule = {}
            for row in csv.DictReader(open(rules_rall_path), delimiter='\t'):
                #NOTE: as of now all the rules are generated using MNX
                #but it may be that other db are used, we are handling this case
                #WARNING: can have multiple products so need to seperate them
                products = {}
                for i in row['Product_IDs'].split('.'):
                    mnxm = self._checkMNXMdeprecated(i)
                    if not mnxm in products:
                        products[mnxm] = 1
                    else:
                        products[mnxm] += 1
                try:
                    #WARNING: one reaction rule can have multiple reactions associated with them
                    #To change when you can set subpaths from the mutliple numbers of
                    #we assume that the reaction rule has multiple unique reactions associated
                    if row['# Rule_ID'] not in rule:
                        rule[row['# Rule_ID']] = {}
                    if row['# Rule_ID'] in rule[row['# Rule_ID']]:
                        self.logger.warning('There is already reaction '+str(row['# Rule_ID'])+' in reaction rule '+str(row['# Rule_ID']))
                    rule[row['# Rule_ID']][row['Reaction_ID']] = {'rule_id': row['# Rule_ID'], 'rule_score': float(row['Score_normalized']), 'reac_id': self._checkMNXRdeprecated(row['Reaction_ID']), 'subs_id': self._checkMNXMdeprecated(row['Substrate_ID']), 'rel_direction': int(row['Rule_relative_direction']), 'left': {self._checkMNXMdeprecated(row['Substrate_ID']): 1}, 'right': products}
                except ValueError:
                    self.logger.error('Problem converting rel_direction: '+str(row['Rule_relative_direction']))
                    self.logger.error('Problem converting rule_score: '+str(row['Score_normalized']))
        except FileNotFoundError as e:
                self.logger.error('Could not read the rules_rall file ('+str(rules_rall_path)+')')
                return {}
        return rule


    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of all cross references for a given chemical id (MNX) to other database id's
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def mnx_chemXref(self, chem_xref_path):
        chemXref = {}
        with open(chem_xref_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = self._checkMNXMdeprecated(row[1])
                    if len(row[0].split(':'))==1:
                        dbName = 'mnx'
                        dbId = row[0]
                    else:
                        dbName = row[0].split(':')[0]
                        dbId = ''.join(row[0].split(':')[1:])
                        if dbName=='deprecated':
                            dbName = 'mnx'
                    #mnx
                    if not mnx in chemXref:
                        chemXref[mnx] = {}
                    if not dbName in chemXref[mnx]:
                        chemXref[mnx][dbName] = []
                    if not dbId in chemXref[mnx][dbName]:
                        chemXref[mnx][dbName].append(dbId)
                    ### DB ###
                    if not dbName in chemXref:
                        chemXref[dbName] = {}
                    if not dbId in chemXref[dbName]:
                        chemXref[dbName][dbId] = mnx
        return chemXref


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
