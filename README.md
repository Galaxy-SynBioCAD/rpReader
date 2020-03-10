# rpReader 

* Docker image: [brsynth/rpreader-standalone](https://hub.docker.com/r/brsynth/rpreader-standalone)

RetroPath2.0 and rp2paths to SBML conversion tool. Takes for input the output of both RetroPath2.0 and rp2paths or RetroPath3.0 and generates a series of enriches SBML files with all the files contained within them. 

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

* Base Docker Image: [brsynth/rpBase](https://hub.docker.com/r/brsynth/rpbase)
* Cache Docker Image: [brsynth/rpCache](https://hub.docker.com/r/brsynth/rpcache)

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

Version 0.1

## Authors

* **Melchior du Lac**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson

## How to cite rpReader?
