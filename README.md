# rp2Reader REST service

REST service to parse the output of RetroPath2.0 and RP2paths to output as a collection of SBML files

## Installing

```
docker build -t brsynth/rpreader-standalone:dev -f Dockerfile .
```

To test untar the test.tar.xz file and run the following command:

```
python run.py -rp2paths_compounds test/rp2paths_compounds.csv -rp2_pathways test/rp2_pathways.csv -rp2paths_pathways test/rp2paths_pathways.csv -outputTar test/test_rpReader.tar
```

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

TODO

## Authors

* **Melchior du Lac**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson
