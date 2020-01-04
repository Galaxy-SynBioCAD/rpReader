#!/bin/bash

python tool_rp2Reader.py -inputTar test_input.tar -outputTar test_output.tar -pathway_id rp_pathway -compartment_id MNXC3
mv test_output.tar results/
