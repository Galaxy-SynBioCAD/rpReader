# rp2Reader REST service 

REST service to parse the output of RetroPath2.0 and RP2paths to output as a collection of SBML files

## Installing

```
docker build -t brsynth/rpfba -f Dockerfile .
```

Run the service

```
docker run --network host -p 8994:8994 brsynth/rpfba
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

