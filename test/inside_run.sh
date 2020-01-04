#!/bin/bash

python tool_rp2Reader.py -outputTar test_output.tar -pathway_id rp_pathway -compartment_id MNXC3 -rp2paths_compounds test_rp2paths_compounds.tsv -rp2paths_pathways test_rp2paths_pathways.csv -rp2_pathways test_rp2_pathways.csv -maxRuleIds 2
mv test_output.tar results/
