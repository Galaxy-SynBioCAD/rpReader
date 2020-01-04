#!/bin/sh

rm -f test_output.tar
docker run -d -p 8888:8888 --name test_rpReader brsynth/rpreader
sleep 10
python tool_rp2Reader.py -outputTar test_output.tar -pathway_id rp_pathway -compartment_id MNXC3 -rp2paths_compounds test_rp2paths_compounds.tsv -rp2paths_pathways test_rp2paths_pathways.csv -rp2_pathways test_rp2_pathways.csv -maxRuleIds 2 -server_url http://0.0.0.0:8888/REST
docker kill test_rpReader
docker rm test_rpReader
