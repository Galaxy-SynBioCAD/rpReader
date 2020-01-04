#!/bin/bash

docker run -v ${PWD}/inside_run.sh:/home/inside_run.sh -v ${PWD}/tool_rpCofactors.py:/home/tool_rp2Reader.py -v ${PWD}/test_rp2_pathways.csv:/home/test_rp2_pathways.csv -v ${PWD}/test_rp2paths_compounds.tsv:/home/test_rp2paths_compounds.tsv -v ${PWD}/test_rp2paths_pathways.csv:/home/test_rp2paths_pathways.csv  -v ${PWD}/results/:/home/results/ --rm brsynth/rpreader /bin/sh /home/inside_run.sh

cp results/test_output.tar .
