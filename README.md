# ColocPy

  ColocPy is  a novel automatic tool for object based colocalization which comes in a simple GUI that requires no programming experience hence rendering it ideal for researchers at all levels. It allows the user to select the images of interest and perform i) an automated and robust image segmentation for region of interest identification (ROI), ii) contemporary pre-processing manipulations and iii) an in depth analysis of the co-localization between objects in two different channels.  
   Outputs include data visualization, to decide upon appropriate analysis and an excel file with all the quantitative analysis. ColocPy is meant to simplify colocalization analysis for experimentalists, however it can also be simply be used to identify cells and particles in each channel idependently.


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Local Installion 

A step by step series of examples that tell you how to get a development env running


```
pip install 'package'
```


### Prerequisites
 To use this software  you need to have Python >= 3.8 and all this additional packages:

* mxnet_mkl
* pyqtgraph
* PyQt5
* numpy (>=1.20.0)
* numba
* scipy
* natsort
* readlif
* matplotlib
* pandas
* read_lif
* cellpose
* scikit-image
* trackpy
* seaborn


Linux and Windows are supported for running this code. For running the graphical user interface you will need Windows. At least 8GB of RAM is required to run this softaware. This softaware has been heavily tested on Windows 10 and Ubuntu 20.04.3. Internet connection is required for installation. Please open an issue if you face any problems.



## Running ColocPy
In order to start the CLI version of ColoPy, open a command line terminal, anaconda prompt is advised. and then type the following:

```
 python ColocPy_CLI.py -i 020721_ch2DDX6.lif -o results -m cyto -f 2 -u 3 --pc1 90 --pc2 80 -t 5 -s 2 -w 1 -b 2 -g 3 -d 2 
```
### Parameter setting

The user should specify the parameters before running. In case the user doesn't specify all of the parameters the default ones will be used. The default parameters have been set as following: tracking space: 5px , sigma: 2 , percentiles(both) : 90 , ch0: 1, ch1: 2, ch2: 3, noise level: 3 and the output files will be called 'result.xlsx/result.pdf' 

| Options | Description |
| ---| :--------:|
| -- i | input file |
| -- o | output file |
| -- m | model: cyto/nuclei |
| -- f | from image |
| -- u  | until image|
| -- pc1 |1st colocalized channel percentile |
| -- pc2 | 2nd colocalized channel percentile |
| -- t | tracking space |
| -- s  | gaussian sigma |
| -- ch0 | nucleus channel |
| -- ch1 | 1st colocalized channel |
| -- ch2 | 2nd colocalized channel |
| -- nl  |WBNS noise level |
	

## Using the GUI

The available GUI version provides more instructions and a page in the documentation here. A brief a summary of the functions provided by the GUI is depicted here:

![alt text](https://github.com/alexiaales/ll/blob/main/format2.PNG "image")

## Comparison of GUI and CLI 

| Parameter| GUI | CLI|
| ---| :--------:| :--------:|
|Input file | yes | yes |
|Output file | no | yes |
|Model | yes | yes |
|Image selection | yes | yes |  
|Percentile def | yes | yes |
|Tracking space | no | yes |
|Gaussian sigma | yes | yes |
|Channel selection | yes | yes |
|WBNS noise level | no | yes |


## Built With

* [Cellpose](https://github.com/MouseLand/cellpose) - A generalist algorithm for cellular segmentation
* [WBNS](https://github.com/NienhausLabKIT/HuepfelM/tree/master/WBNS/python_script) -  Wavelet-based Background and Noise Subtraction
* [Trackpy](https://github.com/soft-matter/trackpy) - Python particle tracking toolkit

## Common issues
* Images must be in RGB and lif format
* Images stacks must not overpass 30 stacks 
* Avoid providing image with great heterogenosity of cells since it might lead to non represenattive results

## Authors

* **Vivian Kalamara**  - [PhD candidate](https://github.com/PurpleBooth)
* **Pantelis Topalis** - [Researcher](https://github.com/PurpleBooth)
* **George Garinis**  - [Professor](https://github.com/PurpleBooth)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


