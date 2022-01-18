# ColocPy

  ColocPy is  a novel automatic tool for object based colocalization which comes in a simple GUI that requires no programming experience hence rendering it ideal for researchers at all levels. It allows the user to select the images of interest and perform i) an automated and robust image segmentation for region of interest identification (ROI), ii) contemporary pre-processing manipulations and iii) an in depth analysis of the co-localization between objects in two different channels.  
	


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Local Installion 

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

Linux and Windows are supported for running this code. For running the graphical user interface you will need Windows. At least 8GB of RAM is required to run thsi softaware. This softaware has been heavily tested on Windows 10 and Ubuntu 20.04.3  Please open an issue if you face any problems.
## Running ColocPy
In order to start the CLI version of ColoPy, open a command line terminal, anaconda prompt is advised. and then type the following:

```
 python ColocPy_CLI.py -i 020721_ch2DDX6.lif -o results -m cyto -f 2 -u 3 --pc1 90 --pc2 80 -t 5 -s 2 -w 1 -b 2 -g 3 -d 2 
```
### Parameter setting

The user should define these parameters before running. In case the 

* -- i < input file>
* -- o < output file>
* -- m < model: cyto/nuclei>
* -- f < from image>
* -- u < until image>
* -- pc1 < 1st colocalized channel percentile >
* -- pc2 < 2nd colocalized channel percentile >
* -- t < tracking space>
* -- s  < gaussian sigma>
* -- cho <nucleus channel>
* -- ch1 <1st colocalized channel>
* -- ch2 <2nd colocalized channel>
* -- nl <WBNS noise level>
	



## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Vivian Kalamara** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)
* **Pantelis Topalis** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)
* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
