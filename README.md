<<<<<<< HEAD
# rpReader 
=======
# rpReader docker 
>>>>>>> standalone-dev-mel

* Docker image: [brsynth/rpreader-standalone](https://hub.docker.com/r/brsynth/rpreader-standalone)

RetroPath2.0 and rp2paths to SBML conversion tool. Takes for input the output of both RetroPath2.0 and rp2paths or RetroPath3.0 and generates a series of enriches SBML files with all the files contained within them. 

<<<<<<< HEAD
=======
## Input

Required:
* **-rp2paths_pathways**: (string) Path to the rp2paths pathways file
* **-rp2paths_compounds**: (string) Path to the rp2paths compounds file
* **-rp2_pathways**: (string) Path to the RetroPath2.0 pathways file

Advanced options:
* **-upper_flux_bound**: (integer, default=9999) Upper flux bound value
* **-lower_flux_bound**: (integer, default=0) Lower flux bound value
* **-maxRuleIds**: (integer, default=2) Number of subpaths per paths
* **-pathway_id**: (string, default=rp_pathway) ID of the heterologous pathway
* **-compartment_id**: (string, default=MNXC3 (i.e. cytoplasm)) Heterologous pathway compartment ID
* **-species_group_id**: (string, default=central_species) ID of the central species, i.e. not cofactors, in the heterologous reactions

## Output

* **-outputTar**: (string) Path to the output tar.xz file

>>>>>>> standalone-dev-mel
## Building the docker

To build the docker locally, run the following command in the root folder of the project:

```
docker build -t brsynth/rpreader-standalone:dev -f Dockerfile .
```

### Running the test

To test untar the test.tar.xz file and run the following command:

```
python run.py -rp2paths_compounds test/rp2paths_compounds.csv -rp2_pathways test/rp2_pathways.csv -rp2paths_pathways test/rp2paths_pathways.csv -outputTar test/test_rpReader.tar
```

## Dependencies

<<<<<<< HEAD
* Base Docker Image: [brsynth/rpBase](https://hub.docker.com/r/brsynth/rpbase)
* Cache Docker Image: [brsynth/rpCache](https://hub.docker.com/r/brsynth/rpcache)
=======
* Base docker image: [brsynth/rpBase](https://hub.docker.com/r/brsynth/rpbase)
* Cache docker image: [brsynth/rpCache](https://hub.docker.com/r/brsynth/rpcache)
>>>>>>> standalone-dev-mel

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Version

<<<<<<< HEAD
Version 0.1
=======
v0.1
>>>>>>> standalone-dev-mel

## Authors

* **Melchior du Lac**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson

## How to cite rpReader?
