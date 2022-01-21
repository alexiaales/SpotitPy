# ColocPy
We present a novel automatic tool for object based colocalization which comes in a simple graphical user interface that requires no programming experience hence rendering it ideal for researchers at all levels. ColocPy allows the user to select the images of interest and perform: i) an automated and robust segmentation of the region of interest (ROI), ii) contemporary pre-processing manipulations for de-noising  and iii) an in depth analysis of the co-localization between objects in two different channels.  

## Introduction

We have developed a free and open source tool implemented in python which allows an object based colocalization analysis with quantifications. The algorithm per se is implemented in both a user friendly plugin and in a command line version which enables more flexible processing options. This tool is applicable to a large variety of biological objects and provides a novel automated approach to extract quantifiable results of object based colocalization analysis while incorporating algorithms for background and noise suppression. The software generates two complementary outputs, one with the quantifications of the selected image and another for visual inspection of the analysis. Together we provide researchers with a user friendly tool that allows to perform a fast, reproducible and robust way of quantifying co-localized particles within a limited time. To illustrates the performance of this software we have applied it on the dataset. 

### Local Installion 

Steps required in order to get ColocPy running on your PC.


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
* Images stacks must not exceed 30 stacks 
* Avoid providing image with great heterogenosity of cells since it might lead to non represenattive results

## Authors

* **Vivian Kalamara**  - [PhD candidate](https://www.researchgate.net/profile/Vivian-Kalamara)
* **George Garinis**  - [Professor](https://www.researchgate.net/profile/George-Garinis)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


