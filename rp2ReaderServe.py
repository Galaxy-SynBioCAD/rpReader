import os
import uuid
import shutil
import json
from datetime import datetime
from flask import Flask, request, jsonify, send_file, abort
from flask_restful import Resource, Api
#from rpviz.main import run
import sys
import io
import tarfile
import libsbml

sys.path.insert(0, '/home/')
import rpReader
import rpCache

#from threading import Thread
from uwsgidecorators import thread
#from tasks import threaded_task, uwsgi_task, spool_task, uwsgi_tasks_task

##############################################
################### REST #####################
##############################################


app = Flask(__name__)
api = Api(app)
#dataFolder = os.path.join( os.path.dirname(__file__),  'data' )

#TODO: test that it works well
#declare the rpReader globally to avoid reading the pickle at every instance

#TODO: test passing the parameters directly
#rpreader = rpReader.rpReader()
rpcache = rpCache.rpCache()
rpcache._loadCache()

def stamp(data, status=1):
    appinfo = {'app': 'rpReader', 'version': '1.0', 
               'author': 'Melchior du Lac',
               'organization': 'BRS',
               'time': datetime.now().isoformat(), 
               'status': status}
    out = appinfo.copy()
    out['data'] = data
    return out


class RestApp(Resource):
    """ REST App."""
    def post(self):
        return jsonify(stamp(None))
    def get(self):
        return jsonify(stamp(None))


class RestQuery(Resource):
    """ REST interface that generates the Design.
        Avoid returning numpy or pandas object in
        order to keep the client lighter.
    """
    def post(self):       
        rp2paths_compounds = request.files['rp2paths_compounds'].read()
        rp2_scope = request.files['rp2_scope'].read()
        rp2paths_outPaths = request.files['rp2paths_outPaths'].read()
        params = json.load(request.files['data'])
        #pass the files to the rpReader
        rpreader = rpReader.rpReader()
        rpreader.deprecatedMNXM_mnxm = rpcache.deprecatedMNXM_mnxm
        rpreader.deprecatedMNXR_mnxr = rpcache.deprecatedMNXR_mnxr
        rpreader.mnxm_strc = rpcache.mnxm_strc
        rpreader.inchikey_mnxm = rpcache.inchikey_mnxm
        rpreader.rr_reactions = rpcache.rr_reactions
        rpreader.chemXref = rpcache.chemXref
        rpreader.compXref = rpcache.compXref
        rpreader.nameCompXref = rpcache.nameCompXref
        rpsbml_paths = rpreader.rp2ToSBML(rp2paths_compounds,
                                      rp2_scope,
                                      rp2paths_outPaths,
                                      None,
                                      int(params['maxRuleIds']),
                                      params['path_id'],
                                      params['compartment_id'])
        #pass the SBML results to a tar
        if rpsbml_paths=={}:
            flask.abort(204)
        outputTar = io.BytesIO()
        tf = tarfile.open(fileobj=outputTar, mode='w:xz')
        for rpsbml_name in rpsbml_paths:
            data = libsbml.writeSBMLToString(rpsbml_paths[rpsbml_name].document).encode('utf-8')
            fiOut = io.BytesIO(data)
            info = tarfile.TarInfo(name=rpsbml_name)
            info.size = len(data)
            tf.addfile(tarinfo=info, fileobj=fiOut)
        ########IMPORTANT######
        tf.close()
        outputTar.seek(0)
        #######################
        return send_file(outputTar, as_attachment=True, attachment_filename='rpReader.tar', mimetype='application/x-tar')


@thread
def uwsgi_rpreader(rp2paths_compounds,
        rp2_scope,
        rp2paths_outPaths,
        tmpOutputFolder,
        maxRuleIds,
        path_id,
        compartment_id):
    return rpsbml_paths


api.add_resource(RestApp, '/REST')
api.add_resource(RestQuery, '/REST/Query')


if __name__== "__main__":
    #debug = os.getenv('USER') == 'mdulac'
    app.run(host="0.0.0.0", port=8997, debug=True, threaded=True)
