# rpReader

* Docker image: [brsynth/rpreader-rest](https://hub.docker.com/r/brsynth/rpreader-rest)

RetroPath2.0 and rp2paths to SBML conversion tool. Takes for input the output of both RetroPath2.0 and rp2paths or RetroPath3.0 and generates a series of enriches SBML files with all the files contained within them. 

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
* **server_url**: (string, default=http://0.0.0.0:8888/REST) IP address of the REST service

## Output

* **-outputTar**: (string) Path to the output tar.xz file

## Dependencies

* Base Docker Image: [brsynth/rpBase](https://hub.docker.com/r/brsynth/rpbase)
* Cache Docker Image: [brsynth/rpCache](https://hub.docker.com/r/brsynth/rpcache)

## Installing

To build the image, use the Dockerfile using the following command:

```
docker build -t brsynth/rpreader-rest -f Dockerfile .
```

To run the service on the localhost as the Galaxy instance run:

```
docker run -p 8888:8888 brsynth/rpreader-rest
```

### Testing

To run the test, untar the test.tar.xz folder and run the following command:

```
python tool_rp2Reader.py -rp2paths_compounds test/test_rp2paths_compounds.tsv -rp2_pathways test/test_rp2_pathways.csv -rp2paths_pathways test/test_rp2paths_pathways.csv -output test/test_rpReader.tar
```

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Version

v0.1

## Authors

* **Melchior du Lac**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson

## How to cite rpReader?
