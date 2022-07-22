# -*- coding: utf-8 -*-


from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter.filedialog import askopenfilename,asksaveasfilename
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
from readlif.reader import LifFile
import pylab
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import seaborn as sns
import PIL
from PIL import Image
import read_lif
import skimage # Library for image manipulation
from skimage import io
from skimage.io import imread # sublibrary from skimage
import trackpy as tp # Library for particle tracking
from scipy.ndimage import gaussian_filter
from collections import Counter, OrderedDict
from pywt import wavedecn, waverecn 
import tifffile as tif  
from scipy.ndimage import gaussian_filter  
from joblib import Parallel, delayed 
import matplotlib.backends.backend_pdf
from matplotlib.backends.backend_pdf import PdfPages
import PyPDF2
from fpdf import FPDF
from PyPDF2 import PdfFileMerger
from reportlab.pdfgen import canvas
from datetime import datetime
import multiprocessing
import seaborn as sns
import tkinter as tk
from cellpose import models
from cellpose import plot
import sys, getopt  

image=[]
number=[]   
sizes=[]
particle_sizes=[]
image=[]
number=[]
fp2=[]
number_of_parti=[]
number_of_part_in_Green=[]
number_of_part_in_Blue=[]
inpu_cell_size=[]
inpu_part_size=[]
negative_cells=[]
cells=[]
clean_div=[]
use_GPU = models.use_gpu()
mode = models.Cellpose(gpu=use_GPU, model_type='cyto') 
pdf_files=[]


def pdf_output(img):    
    pdfs = ['GFG.pdf', '1st Channel stacked.pdf', '2nd Channel stacked.pdf', '3rd Channel stacked.pdf', 'Masked nucleus.pdf','Masked cytoplasm.pdf','Calculation area.pdf','Filtered channel 2 image.pdf','Filtered channel 3 image.pdf','Identified particles in 2nd channel.pdf','Identified particles in 3rd channel.pdf']
    merger = PdfFileMerger()    
    for pdf in pdfs:
        merger.append(pdf)
    merger.write(str(img)+'.pdf')
    pdfs_names=str(img)+'.pdf'
    pdf_files.append(str(pdfs_names))
    merger.close()
    pdfs = ['GFG.pdf', '1st Channel stacked.pdf', '2nd Channel stacked.pdf', '3rd Channel stacked.pdf','Masked nucleus.pdf','Masked cytoplasm.pdf','Calculation area.pdf','Filtered channel 2 image.pdf','Filtered channel 3 image.pdf','Identified particles in 2nd channel.pdf','Identified particles in 3rd channel.pdf']
    for f in pdfs:
        os.remove(os.path.join('./', f))
        
def pdf_output2(img):    
    pdfs = ['GFG.pdf', '1st Channel stacked.pdf', '2nd Channel stacked.pdf', '3rd Channel stacked.pdf', 'Masked nucleus.pdf','Calculation area.pdf','Filtered channel 2 image.pdf','Filtered channel 3 image.pdf','Identified particles in 2nd channel.pdf','Identified particles in 3rd channel.pdf']
    merger = PdfFileMerger()    
    for pdf in pdfs:
        merger.append(pdf)
    merger.write(str(img)+'.pdf')
    pdfs_names=str(img)+'.pdf'
    pdf_files.append(str(pdfs_names))
    merger.close()
    pdfs = ['GFG.pdf', '1st Channel stacked.pdf', '2nd Channel stacked.pdf', '3rd Channel stacked.pdf','Masked nucleus.pdf','Masked cytoplasm.pdf','Calculation area.pdf','Filtered channel 2 image.pdf','Filtered channel 3 image.pdf','Identified particles in 2nd channel.pdf','Identified particles in 3rd channel.pdf']
    for f in pdfs:
        os.remove(os.path.join('./', f))
    
def final_merging(counter):
    merger = PdfFileMerger()
    dateTimeObj = datetime.now()
    date='This file was created the: '+str(dateTimeObj.year)+'/'+str(dateTimeObj.month)+'/'+str(dateTimeObj.day)
    title='IMAGE ANALYSIS REPORT'
    per='Persentile: '+str(persentile)
    filt='Sigma Gaussian filter: '+str(sigm)
    num='Number of analyzed images: '+str(counter)
    c = canvas.Canvas("Intro.pdf")
    c.setFont("Courier", 9) #choose your font type and font size
    c.drawString(390, 820, date)
    c.setFont("Helvetica", 26) #choose your font type and font size
    c.drawString(130, 620, title)
    c.setFont("Courier", 9) #choose your font type and font size
    c.drawString(30, 80, per)
    c.drawString(30, 50, filt)
    c.save()  
  #  os.remove(os.path.join('./', 'Saturated image'))
  #  os.remove(os.path.join('./', 'Channel 2 Image'))
  #  os.remove(os.path.join('./', 'Channel 3 Image'))
  #  pdf_files.append(str('Intro.pdf'))
    for pdf in pdf_files:
        merger.append(pdf)
    merger.write(str(outputfile)+'.pdf')
    merger.close()
    for x in pdf_files:
        os.remove(os.path.join('./', x))
    print('Analysis completed')

def plot(image,color,plot_title):     
    plt.figure(figsize=(10,10))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)
    plt.savefig(plot_title+'.pdf') 
    
def ploting(image,color,plot_title):     
    plt.figure(figsize=(10,10))
    plt.imshow(image,cmap=color)
    plt.savefig(plot_title+'.pdf') 
    
def plotting(image,color,plot_title): 
    plt.figure(figsize=(10,10))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)
    plt.savefig(plot_title+'.pdf') 
                
def plotting2(image,color,plot_title): 
    plt.figure(figsize=(10,10))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)
    plt.savefig(plot_title+'.pdf')       

def plotting3(image,color,plot_title):
    plt.figure(figsize=(10,10))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)
    plt.savefig(plot_title+'.pdf') 

def plotting4(image,color,plot_title):     
    plt.figure(figsize=(10,10))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)
    plt.savefig(plot_title+'.pdf') 
    
def plotting5(image,color,plot_title):     
    plt.figure(figsize=(10,10))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)
    plt.savefig(plot_title+'.pdf') 
 
def stacking_step(z_list_0,z_list_1,z_list_2):
    z_st_0=[]
    for x in range(len(z_list_0)):
        pic=np.array(z_list_0[x])
        z_st_0.append(pic)   
    if len(z_st_0)==19:
        global z_0
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13],z_st_0[14],z_st_0[15],z_st_0[16],z_st_0[17],z_st_0[18]])              
    elif len(z_st_0)==30:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13],z_st_0[14],z_st_0[15],z_st_0[16],z_st_0[17],z_st_0[18],z_st_0[19],z_st_0[20],z_st_0[21],z_st_0[22],z_st_0[23],z_st_0[24],z_st_0[25],z_st_0[26],z_st_0[27],z_st_0[28],z_st_0[29]])
    elif len(z_st_0)==29:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13],z_st_0[14],z_st_0[15],z_st_0[16],z_st_0[17],z_st_0[18],z_st_0[19],z_st_0[20],z_st_0[21],z_st_0[22],z_st_0[23],z_st_0[24],z_st_0[25],z_st_0[26],z_st_0[27],z_st_0[28]])
    elif len(z_st_0)==28:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13],z_st_0[14],z_st_0[15],z_st_0[16],z_st_0[17],z_st_0[18],z_st_0[19],z_st_0[20],z_st_0[21],z_st_0[22],z_st_0[23],z_st_0[24],z_st_0[25],z_st_0[26],z_st_0[27]])
    elif len(z_st_0)==27:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13],z_st_0[14],z_st_0[15],z_st_0[16],z_st_0[17],z_st_0[18],z_st_0[19],z_st_0[20],z_st_0[21],z_st_0[22],z_st_0[23],z_st_0[24],z_st_0[25],z_st_0[26]])
    elif len(z_st_0)==26:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13],z_st_0[14],z_st_0[15],z_st_0[16],z_st_0[17],z_st_0[18],z_st_0[19],z_st_0[20],z_st_0[21],z_st_0[22],z_st_0[23],z_st_0[24],z_st_0[25]])
    elif len(z_st_0)==25:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13],z_st_0[14],z_st_0[15],z_st_0[16],z_st_0[17],z_st_0[18],z_st_0[19],z_st_0[20],z_st_0[21],z_st_0[22],z_st_0[23],z_st_0[24]])
    elif len(z_st_0)==24:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13],z_st_0[14],z_st_0[15],z_st_0[16],z_st_0[17],z_st_0[18],z_st_0[19],z_st_0[20],z_st_0[21],z_st_0[22],z_st_0[23]])
    elif len(z_st_0)==23:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13],z_st_0[14],z_st_0[15],z_st_0[16],z_st_0[17],z_st_0[18],z_st_0[19],z_st_0[20],z_st_0[21],z_st_0[22]])
    elif len(z_st_0)==22:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13],z_st_0[14],z_st_0[15],z_st_0[16],z_st_0[17],z_st_0[18],z_st_0[19],z_st_0[20],z_st_0[21]])
    elif len(z_st_0)==21:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13],z_st_0[14],z_st_0[15],z_st_0[16],z_st_0[17],z_st_0[18],z_st_0[19],z_st_0[20]])
    elif len(z_st_0)==20:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13],z_st_0[14],z_st_0[15],z_st_0[16],z_st_0[17],z_st_0[18],z_st_0[19]])    
    elif len(z_st_0)==18:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13],z_st_0[14],z_st_0[15],z_st_0[16],z_st_0[17]])              
    elif len(z_st_0)==17:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13],z_st_0[14],z_st_0[15],z_st_0[16]])          
    elif len(z_st_0)==16:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13],z_st_0[14],z_st_0[15]])      
    elif len(z_st_0)==15:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13],z_st_0[14]])      
    elif len(z_st_0)==14:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12],z_st_0[13]])  
    elif len(z_st_0)==13:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11],z_st_0[12]])  
    elif len(z_st_0)==12:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10],z_st_0[11]])  
    elif len(z_st_0)==11:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9],z_st_0[10]]) 
    elif len(z_st_0)==10:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8],z_st_0[9]]) 
    elif len(z_st_0)==9:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7],z_st_0[8]]) 
    elif len(z_st_0)==8:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6],z_st_0[7]]) 
    elif len(z_st_0)==7:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5],z_st_0[6]]) 
    elif len(z_st_0)==6:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4],z_st_0[5]]) 
    elif len(z_st_0)==5:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3],z_st_0[4]]) 
    elif len(z_st_0)==4:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2],z_st_0[3]]) 
    elif len(z_st_0)==3:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1],z_st_0[2]]) 
    elif len(z_st_0)==2:
        z_0=np.maximum.reduce([z_st_0[0],z_st_0[1]])  
    elif len(z_st_0)==1:
        z_0=np.maximum.reduce([z_st_0[0]])
        
    plot(z_0,'Reds','1st Channel stacked')
        
    z_st_1=[]
    for x in range(len(z_list_1)):
        pic=np.array(z_list_1[x])
        z_st_1.append(pic)   
 
                
    if len(z_st_1)==19:
        global z_1
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13],z_st_1[14],z_st_1[15],z_st_1[16],z_st_1[17],z_st_1[18]])              
    elif len(z_st_1)==30:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13],z_st_1[14],z_st_1[15],z_st_1[16],z_st_1[17],z_st_1[18],z_st_1[19],z_st_1[20],z_st_1[21],z_st_1[22],z_st_1[23],z_st_1[24],z_st_1[25],z_st_1[26],z_st_1[27],z_st_1[28],z_st_1[29]])
    elif len(z_st_1)==29:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13],z_st_1[14],z_st_1[15],z_st_1[16],z_st_1[17],z_st_1[18],z_st_1[19],z_st_1[20],z_st_1[21],z_st_1[22],z_st_1[23],z_st_1[24],z_st_1[25],z_st_1[26],z_st_1[27],z_st_1[28]])
    elif len(z_st_1)==28:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13],z_st_1[14],z_st_1[15],z_st_1[16],z_st_1[17],z_st_1[18],z_st_1[19],z_st_1[20],z_st_1[21],z_st_1[22],z_st_1[23],z_st_1[24],z_st_1[25],z_st_1[26],z_st_1[27]])
    elif len(z_st_1)==27:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13],z_st_1[14],z_st_1[15],z_st_1[16],z_st_1[17],z_st_1[18],z_st_1[19],z_st_1[20],z_st_1[21],z_st_1[22],z_st_1[23],z_st_1[24],z_st_1[25],z_st_1[26]])
    elif len(z_st_1)==26:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13],z_st_1[14],z_st_1[15],z_st_1[16],z_st_1[17],z_st_1[18],z_st_1[19],z_st_1[20],z_st_1[21],z_st_1[22],z_st_1[23],z_st_1[24],z_st_1[25]])
    elif len(z_st_1)==25:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13],z_st_1[14],z_st_1[15],z_st_1[16],z_st_1[17],z_st_1[18],z_st_1[19],z_st_1[20],z_st_1[21],z_st_1[22],z_st_1[23],z_st_1[24]])
    elif len(z_st_1)==24:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13],z_st_1[14],z_st_1[15],z_st_1[16],z_st_1[17],z_st_1[18],z_st_1[19],z_st_1[20],z_st_1[21],z_st_1[22],z_st_1[23]])
    elif len(z_st_1)==23:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13],z_st_1[14],z_st_1[15],z_st_1[16],z_st_1[17],z_st_1[18],z_st_1[19],z_st_1[20],z_st_1[21],z_st_1[22]])
    elif len(z_st_1)==22:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13],z_st_1[14],z_st_1[15],z_st_1[16],z_st_1[17],z_st_1[18],z_st_1[19],z_st_1[20],z_st_1[21]])
    elif len(z_st_1)==21:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13],z_st_1[14],z_st_1[15],z_st_1[16],z_st_1[17],z_st_1[18],z_st_1[19],z_st_1[20]])
    elif len(z_st_1)==20:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13],z_st_1[14],z_st_1[15],z_st_1[16],z_st_1[17],z_st_1[18],z_st_1[19]])        
    elif len(z_st_1)==18:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13],z_st_1[14],z_st_1[15],z_st_1[16],z_st_1[17]])              
    elif len(z_st_1)==17:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13],z_st_1[14],z_st_1[15],z_st_1[16]])          
    elif len(z_st_1)==16:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13],z_st_1[14],z_st_1[15]])      
    elif len(z_st_1)==15:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13],z_st_1[14]])      
    elif len(z_st_1)==14:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12],z_st_1[13]])  
    elif len(z_st_1)==13:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11],z_st_1[12]])  
    elif len(z_st_1)==12:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10],z_st_1[11]])  
    elif len(z_st_1)==11:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9],z_st_1[10]]) 
    elif len(z_st_1)==10:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8],z_st_1[9]]) 
    elif len(z_st_1)==9:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7],z_st_1[8]]) 
    elif len(z_st_1)==8:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6],z_st_1[7]]) 
    elif len(z_st_1)==7:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5],z_st_1[6]]) 
    elif len(z_st_1)==6:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4],z_st_1[5]]) 
    elif len(z_st_1)==5:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3],z_st_1[4]]) 
    elif len(z_st_1)==4:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2],z_st_1[3]]) 
    elif len(z_st_1)==3:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1],z_st_1[2]]) 
    elif len(z_st_1)==2:
        z_1=np.maximum.reduce([z_st_1[0],z_st_1[1]]) 
    elif len(z_st_1)==1:
        z_1=np.maximum.reduce([z_st_1[0]]) 
        
    plot(z_1,'Greens','2nd Channel stacked')
            
    
    z_st_2=[]
    for x in range(len(z_list_2)):
        pic=np.array(z_list_2[x])
        z_st_2.append(pic)   
            
    if len(z_st_2)==19:
        global z_2
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13],z_st_2[14],z_st_2[15],z_st_2[16],z_st_2[17],z_st_2[18]])          
    elif len(z_st_2)==30:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13],z_st_2[14],z_st_2[15],z_st_2[16],z_st_2[17],z_st_2[18],z_st_2[19],z_st_2[20],z_st_2[21],z_st_2[22],z_st_2[23],z_st_2[24],z_st_2[25],z_st_2[26],z_st_2[27],z_st_2[28],z_st_2[29]])
    elif len(z_st_2)==29:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13],z_st_2[14],z_st_2[15],z_st_2[16],z_st_2[17],z_st_2[18],z_st_2[19],z_st_2[20],z_st_2[21],z_st_2[22],z_st_2[23],z_st_2[24],z_st_2[25],z_st_2[26],z_st_2[27],z_st_2[28]])
    elif len(z_st_2)==28:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13],z_st_2[14],z_st_2[15],z_st_2[16],z_st_2[17],z_st_2[18],z_st_2[19],z_st_2[20],z_st_2[21],z_st_2[22],z_st_2[23],z_st_2[24],z_st_2[25],z_st_2[26],z_st_2[27]])
    elif len(z_st_2)==27:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13],z_st_2[14],z_st_2[15],z_st_2[16],z_st_2[17],z_st_2[18],z_st_2[19],z_st_2[20],z_st_2[21],z_st_2[22],z_st_2[23],z_st_2[24],z_st_2[25],z_st_2[26]])
    elif len(z_st_2)==26:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13],z_st_2[14],z_st_2[15],z_st_2[16],z_st_2[17],z_st_2[18],z_st_2[19],z_st_2[20],z_st_2[21],z_st_2[22],z_st_2[23],z_st_2[24],z_st_2[25]])
    elif len(z_st_2)==25:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13],z_st_2[14],z_st_2[15],z_st_2[16],z_st_2[17],z_st_2[18],z_st_2[19],z_st_2[20],z_st_2[21],z_st_2[22],z_st_2[23],z_st_2[24]])
    elif len(z_st_2)==24:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13],z_st_2[14],z_st_2[15],z_st_2[16],z_st_2[17],z_st_2[18],z_st_2[19],z_st_2[20],z_st_2[21],z_st_2[22],z_st_2[23]])
    elif len(z_st_2)==23:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13],z_st_2[14],z_st_2[15],z_st_2[16],z_st_2[17],z_st_2[18],z_st_2[19],z_st_2[20],z_st_2[21],z_st_2[22]])
    elif len(z_st_2)==22:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13],z_st_2[14],z_st_2[15],z_st_2[16],z_st_2[17],z_st_2[18],z_st_2[19],z_st_2[20],z_st_2[21]])
    elif len(z_st_2)==21:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13],z_st_2[14],z_st_2[15],z_st_2[16],z_st_2[17],z_st_2[18],z_st_2[19],z_st_2[20]])
    elif len(z_st_2)==20:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13],z_st_2[14],z_st_2[15],z_st_2[16],z_st_2[17],z_st_2[18],z_st_2[19]])        
    elif len(z_st_2)==18:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13],z_st_2[14],z_st_2[15],z_st_2[16],z_st_2[17]])              
    elif len(z_st_2)==17:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13],z_st_2[14],z_st_2[15],z_st_2[16]])          
    elif len(z_st_2)==16:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13],z_st_2[14],z_st_2[15]])      
    elif len(z_st_2)==15:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13],z_st_2[14]])      
    elif len(z_st_2)==14:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12],z_st_2[13]])  
    elif len(z_st_2)==13:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11],z_st_2[12]])  
    elif len(z_st_2)==12:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10],z_st_2[11]])  
    elif len(z_st_2)==11:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9],z_st_2[10]]) 
    elif len(z_st_2)==10:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8],z_st_2[9]]) 
    elif len(z_st_2)==9:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7],z_st_2[8]]) 
    elif len(z_st_2)==8:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6],z_st_2[7]]) 
    elif len(z_st_2)==7:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5],z_st_2[6]]) 
    elif len(z_st_2)==6:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4],z_st_2[5]]) 
    elif len(z_st_2)==5:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3],z_st_2[4]]) 
    elif len(z_st_2)==4:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2],z_st_2[3]]) 
    elif len(z_st_2)==3:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1],z_st_2[2]]) 
    elif len(z_st_2)==2:
        z_2=np.maximum.reduce([z_st_2[0],z_st_2[1]]) 
    elif len(z_st_2)==1:
        z_2=np.maximum.reduce([z_st_2[0]]) 
        
    plot(z_2,'Blues','3rd Channel stacked')
        
               
    return z_0,z_1,z_2      
        
def saturated_images(z_list_0,z_list_2):
    #Images with maximum intensity    
    z_st_0=[]
    for x in range(len(z_list_0)):
        pic=np.array(z_list_0[x])
        z_st_0.append(pic)   
    if len(z_st_0)==10:
        global z_3
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]
    elif len(z_st_0)==30:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]+z_st_0[14]+z_st_0[15]+z_st_0[16]+z_st_0[17]+z_st_0[18]+z_st_0[19]+z_st_0[20]+z_st_0[21]+z_st_0[22]+z_st_0[23]+z_st_0[24]+z_st_0[25]+z_st_0[26]+z_st_0[27]+z_st_0[28]+z_st_0[29]
    elif len(z_st_0)==29:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]+z_st_0[14]+z_st_0[15]+z_st_0[16]+z_st_0[17]+z_st_0[18]+z_st_0[19]+z_st_0[20]+z_st_0[21]+z_st_0[22]+z_st_0[23]+z_st_0[24]+z_st_0[25]+z_st_0[26]+z_st_0[27]+z_st_0[28]
    elif len(z_st_0)==28:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]+z_st_0[14]+z_st_0[15]+z_st_0[16]+z_st_0[17]+z_st_0[18]+z_st_0[19]+z_st_0[20]+z_st_0[21]+z_st_0[22]+z_st_0[23]+z_st_0[24]+z_st_0[25]+z_st_0[26]+z_st_0[27]
    elif len(z_st_0)==27:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]+z_st_0[14]+z_st_0[15]+z_st_0[16]+z_st_0[17]+z_st_0[18]+z_st_0[19]+z_st_0[20]+z_st_0[21]+z_st_0[22]+z_st_0[23]+z_st_0[24]+z_st_0[25]+z_st_0[26]
    elif len(z_st_0)==26:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]+z_st_0[14]+z_st_0[15]+z_st_0[16]+z_st_0[17]+z_st_0[18]+z_st_0[19]+z_st_0[20]+z_st_0[21]+z_st_0[22]+z_st_0[23]+z_st_0[24]+z_st_0[25]
    elif len(z_st_0)==25:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]+z_st_0[14]+z_st_0[15]+z_st_0[16]+z_st_0[17]+z_st_0[18]+z_st_0[19]+z_st_0[20]+z_st_0[21]+z_st_0[22]+z_st_0[23]+z_st_0[24]
    elif len(z_st_0)==24:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]+z_st_0[14]+z_st_0[15]+z_st_0[16]+z_st_0[17]+z_st_0[18]+z_st_0[19]+z_st_0[20]+z_st_0[21]+z_st_0[22]+z_st_0[23]
    elif len(z_st_0)==23:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]+z_st_0[14]+z_st_0[15]+z_st_0[16]+z_st_0[17]+z_st_0[18]+z_st_0[19]+z_st_0[20]+z_st_0[21]+z_st_0[22]
    elif len(z_st_0)==22:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]+z_st_0[14]+z_st_0[15]+z_st_0[16]+z_st_0[17]+z_st_0[18]+z_st_0[19]+z_st_0[20]+z_st_0[21]
    elif len(z_st_0)==21:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]+z_st_0[14]+z_st_0[15]+z_st_0[16]+z_st_0[17]+z_st_0[18]+z_st_0[19]+z_st_0[20]
    elif len(z_st_0)==20:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]+z_st_0[14]+z_st_0[15]+z_st_0[16]+z_st_0[17]+z_st_0[18]+z_st_0[19]    
    elif len(z_st_0)==19:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]+z_st_0[14]+z_st_0[15]+z_st_0[16]+z_st_0[17]+z_st_0[18]
    elif len(z_st_0)==18:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]+z_st_0[14]+z_st_0[15]+z_st_0[16]+z_st_0[17]
    elif len(z_st_0)==17:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]+z_st_0[14]+z_st_0[15]+z_st_0[16]
    elif len(z_st_0)==16:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]+z_st_0[14]+z_st_0[15]
    elif len(z_st_0)==15:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]+z_st_0[14]
    elif len(z_st_0)==14:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]+z_st_0[13]
    elif len(z_st_0)==13:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]+z_st_0[12]
    elif len(z_st_0)==12:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]+z_st_0[11]   
    elif len(z_st_0)==11:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]+z_st_0[9]+z_st_0[10]  
    elif len(z_st_0)==9:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]+z_st_0[8]
    elif len(z_st_0)==8:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]+z_st_0[7]
    elif len(z_st_0)==7:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]+z_st_0[6]
    elif len(z_st_0)==6:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]+z_st_0[5]
    elif len(z_st_0)==5:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]+z_st_0[4]
    elif len(z_st_0)==4:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]+z_st_0[3]
    elif len(z_st_0)==3:
        z_3=z_st_0[0]+z_st_0[1]+z_st_0[2]
    elif len(z_st_0)==2:
        z_3=z_st_0[0]+z_st_0[1]
    elif len(z_st_0)==1:
        z_3=z_st_0[0]  
    
        
    
    z_st_2=[]
    for x in range(len(z_list_2)):
        pic=np.array(z_list_2[x])
        z_st_2.append(pic)    
           
    if len(z_st_2)==10:
        global z_4
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]   
    elif len(z_st_2)==30:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]+z_st_2[14]+z_st_2[15]+z_st_2[16]+z_st_2[17]+z_st_2[18]+z_st_2[19]+z_st_2[20]+z_st_2[21]+z_st_2[22]+z_st_2[23]+z_st_2[24]+z_st_2[25]+z_st_2[26]+z_st_2[27]+z_st_2[28]+z_st_2[29]
    elif len(z_st_2)==29:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]+z_st_2[14]+z_st_2[15]+z_st_2[16]+z_st_2[17]+z_st_2[18]+z_st_2[19]+z_st_2[20]+z_st_2[21]+z_st_2[22]+z_st_2[23]+z_st_2[24]+z_st_2[25]+z_st_2[26]+z_st_2[27]+z_st_2[28]
    elif len(z_st_2)==28:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]+z_st_2[14]+z_st_2[15]+z_st_2[16]+z_st_2[17]+z_st_2[18]+z_st_2[19]+z_st_2[20]+z_st_2[21]+z_st_2[22]+z_st_2[23]+z_st_2[24]+z_st_2[25]+z_st_2[26]+z_st_2[27]
    elif len(z_st_2)==27:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]+z_st_2[14]+z_st_2[15]+z_st_2[16]+z_st_2[17]+z_st_2[18]+z_st_2[19]+z_st_2[20]+z_st_2[21]+z_st_2[22]+z_st_2[23]+z_st_2[24]+z_st_2[25]+z_st_2[26]
    elif len(z_st_2)==26:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]+z_st_2[14]+z_st_2[15]+z_st_2[16]+z_st_2[17]+z_st_2[18]+z_st_2[19]+z_st_2[20]+z_st_2[21]+z_st_2[22]+z_st_2[23]+z_st_2[24]+z_st_2[25]
    elif len(z_st_2)==25:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]+z_st_2[14]+z_st_2[15]+z_st_2[16]+z_st_2[17]+z_st_2[18]+z_st_2[19]+z_st_2[20]+z_st_2[21]+z_st_2[22]+z_st_2[23]+z_st_2[24]
    elif len(z_st_2)==24:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]+z_st_2[14]+z_st_2[15]+z_st_2[16]+z_st_2[17]+z_st_2[18]+z_st_2[19]+z_st_2[20]+z_st_2[21]+z_st_2[22]+z_st_2[23]
    elif len(z_st_2)==23:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]+z_st_2[14]+z_st_2[15]+z_st_2[16]+z_st_2[17]+z_st_2[18]+z_st_2[19]+z_st_2[20]+z_st_2[21]+z_st_2[22]
    elif len(z_st_2)==22:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]+z_st_2[14]+z_st_2[15]+z_st_2[16]+z_st_2[17]+z_st_2[18]+z_st_2[19]+z_st_2[20]+z_st_2[21]
    elif len(z_st_2)==21:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]+z_st_2[14]+z_st_2[15]+z_st_2[16]+z_st_2[17]+z_st_2[18]+z_st_2[19]+z_st_2[20]
    elif len(z_st_2)==20:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]+z_st_2[14]+z_st_2[15]+z_st_2[16]+z_st_2[17]+z_st_2[18]+z_st_2[19]      
    elif len(z_st_2)==19:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]+z_st_2[14]+z_st_2[15]+z_st_2[16]+z_st_2[17]+z_st_2[18]
    elif len(z_st_2)==18:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]+z_st_2[14]+z_st_2[15]+z_st_2[16]+z_st_2[17]
    elif len(z_st_2)==17:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]+z_st_2[14]+z_st_2[15]+z_st_2[16]
    elif len(z_st_2)==16:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]+z_st_2[14]+z_st_2[15]
    elif len(z_st_2)==15:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]+z_st_2[14]
    elif len(z_st_2)==14:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]+z_st_2[13]
    elif len(z_st_2)==13:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]+z_st_2[12]
    elif len(z_st_2)==12:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]+z_st_2[11]   
    elif len(z_st_2)==11:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]+z_st_2[9]+z_st_2[10]  
    elif len(z_st_2)==9:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]+z_st_2[8]
    elif len(z_st_2)==8:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]+z_st_2[7]
    elif len(z_st_2)==7:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]+z_st_2[6]
    elif len(z_st_2)==6:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]+z_st_2[5]
    elif len(z_st_2)==5:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]+z_st_2[4]
    elif len(z_st_2)==4:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]+z_st_2[3]
    elif len(z_st_2)==3:
        z_4=z_st_2[0]+z_st_2[1]+z_st_2[2]
    elif len(z_st_2)==2:
        z_4=z_st_2[0]+z_st_2[1]
    elif len(z_st_2)==1:
        z_4=z_st_2[0]   
      
    return z_3,z_4 

def wavelet_based_BG_subtraction(image,num_levels,noise_lvl):

  coeffs = wavedecn(image, 'db1', level=None) #decomposition
  coeffs2 = coeffs.copy()
  
  for BGlvl in range(1, num_levels):
      coeffs[-BGlvl] = {k: np.zeros_like(v) for k, v in coeffs[-BGlvl].items()} #set lvl 1 details  to zero
  
  Background = waverecn(coeffs, 'db1') #reconstruction
  del coeffs
  BG_unfiltered = Background
  Background = gaussian_filter(Background, sigma=2**num_levels) #gaussian filter sigma = 2^#lvls 
  
  coeffs2[0] = np.ones_like(coeffs2[0]) #set approx to one (constant)
  for lvl in range(1, np.size(coeffs2)-noise_lvl):
      coeffs2[lvl] = {k: np.zeros_like(v) for k, v in coeffs2[lvl].items()} #keep first detail lvl only
  Noise = waverecn(coeffs2, 'db1') #reconstruction
  del coeffs2
  
  return Background, Noise, BG_unfiltered

def WBNS(input_image):
# Wavelet-based background and noise subtraction for fluorescence microscopy images

    resolution_px = 6
    #insert the noise level. If resolution_px > 6 then noise_lvl = 2 may be better 
    noise_lvl = int(nl) #default = 1
    #number of levels for background estimate
    num_levels = np.uint16(np.ceil(np.log2(resolution_px)))
    #read image file adjust shape if neccessary (padding) and plot
    #padding=extending the area of which a convolutional neural network processes an image.
    image = input_image
    img_type = image.dtype
    image = np.array(image,dtype = 'float32')

    #image = np.array(io.imread(os.path.join(data_dir, file)),dtype = 'float32')
    if np.ndim(image) == 2:
        shape = np.shape(image)
        image = np.reshape(image, [1, shape[0], shape[1]])
    shape = np.shape(image)
    if shape[1] % 2 != 0:
        image = np.pad(image,((0,0), (0,1), (0, 0)), 'edge')
        pad_1 = True
    else:
        pad_1 = False
    if shape[2] % 2 != 0:
         image = np.pad(image,((0,0), (0,0), (0, 1)), 'edge')
         pad_2 = True
    else:
         pad_2 = False
         
    #extract background and noise
    num_cores = multiprocessing.cpu_count() #number of cores on your CPU
    res = Parallel(n_jobs=num_cores,max_nbytes=None)(delayed(wavelet_based_BG_subtraction)(image[slice],num_levels, noise_lvl) for slice in range(np.size(image,0)))
    Background, Noise, BG_unfiltered = zip(*res)
    #convert to float64 numpy array
    Noise = np.asarray(Noise,dtype = 'float32')
    Background = np.asarray(Background,dtype = 'float32')
    BG_unfiltered = np.asarray(BG_unfiltered,dtype = 'float32')
    #undo padding
    if pad_1:
        image = image[:,:-1,:]
        Noise = Noise[:,:-1,:]
        Background = Background[:,:-1,:]
        BG_unfiltered = BG_unfiltered[:,:-1,:]
    if pad_2:
        image = image[:,:,:-1]
        Noise = Noise[:,:,:-1]
        Background = Background[:,:,:-1]
        BG_unfiltered = BG_unfiltered[:,:,:-1]
    
    #save unfiltered BG
    BG_unfiltered = np.asarray(BG_unfiltered,dtype=img_type.name)
    #save and plot filtered BG
    Background = np.asarray(Background,dtype=img_type.name)
    #plt.figure()
    #imgplot = plt.imshow(np.amax(Background,0), cmap='Blues')
    #subtract BG only
    result = image - Background
    result[result<0] = 0 #positivity constraint
    #save and plot noisy signal
    result = np.asarray(result,dtype=img_type.name)
    noisy_sig = result
    #correct noise
    Noise[Noise<0] = 0 #positivity constraint
    noise_threshold = np.mean(Noise)+2*np.std(Noise)
    Noise[Noise>noise_threshold] = noise_threshold #2 sigma threshold reduces artifacts
    #subtract Noise
    result = image - Background
    result = result - Noise
    result[result<0] = 0 #positivity constraint
    #save noise
    Noise = np.asarray(Noise,dtype=img_type.name)
    #save result
    result = np.asarray(result,dtype=img_type.name)
    #Plot result
    #plt.figure(figsize=(10,10))
    #imgplot = plt.imshow(np.amax(result,0), cmap='Blues')
    
    return result

def masking(nucleus,cytoplasm,image):
    import matplotlib.pyplot as plt
    cell_s=sizes[int(image)]
    masks, flows, styles, diams = mode.eval(nucleus, diameter=int(cell_s), flow_threshold=None, channels=[0,0])
    masks= gaussian_filter(masks, sigma=int(sigm))
    
    plotting3(masks,'Reds','Masked nucleus')

    #model = models.Cellpose(gpu=use_GPU, model_type='cyto') 
    masksc, flowsc, stylesc, diamsc = mode.eval(cytoplasm, diameter=int(cell_s)*2, flow_threshold=None, channels=[0,0])
    masksc= gaussian_filter(masksc, sigma=int(sigm))
    plotting4(masksc,'Blues','Masked cytoplasm')
    
    masksc_copy=masksc.copy()
    masksc[masksc>1]=1
    masks_copy = masks.copy()
    masks[masks==0]=266
    masks[masks<266]=0
    
    masking_preview = np.einsum('jk,jk->jk',masksc_copy, masks) # remove the nucleus for  noise 
    plotting4(masking_preview, 'Oranges', 'Calculation area')
    
    return masks,masksc_copy,masks_copy,masksc

def masking2(nucleus,cytoplasm,image):
 
    cell_s=sizes[int(image)] 
    masks, flows, styles, diams = mode.eval(nucleus, diameter=int(cell_s), flow_threshold=None, channels=[0,0])
    masks_copy = masks.copy()
    masks= gaussian_filter(masks, sigma=int(sigm))
    plotting(masks,'Oranges','Masked nucleus')
    
    masksc, flowsc, stylesc, diamsc = mode.eval(cytoplasm, diameter=int(cell_s)*1.5, flow_threshold=None, channels=[0,0])
    masksc= gaussian_filter(masksc, sigma=int(sigm))
    plotting(masksc,'Oranges','Masked cytoplasm')
    
    masks[masks>1]=1
    plotting(masks,'Oranges','Calculation area')
    
    return masks,masksc,masks_copy

def quantitative_analysis(nucl_image):
    #find how many cell identifiers exist in my array and see how many times each identifier occurs
    counts=[]
    for x in range(len(np.unique(nucl_image))):
        if x == 0:
            continue
        count = np.count_nonzero(nucl_image == int(x))
        counts.append(count)
    
    # I find the biggest cell and see how far away the other ones are in oredr to create the exact number of cells so that i can devide in the end    
    cells=0
    if len(counts)>5:
        avg=sum(counts)/len(counts)
        for x in counts:
            percentage=int(x)/int(avg)*100
            if percentage > 50:
                cells+=100
            else:
                cells+=int(percentage)
    else:           
        for x in counts:
            percentage=int(x)/int(max(counts))*100
            if percentage > 49:
                cells+=100
            else:
                cells+=int(percentage)
    
    number_of_cells=cells/100 
    
    
    return(counts,number_of_cells)
   
def tracking3(img1,img2,number_of_cells,mask_cyto,image,img):
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    particle_size=particle_sizes[int(image)]
    ts=particle_sizes[int(image)]
    f1 = tp.locate(img1,int(particle_size),percentile=int((persentile)))
    plt.figure(figsize=(50,50))
    plt.title('Identified particles in 2nd channel')
    f=tp.annotate(f1, img1)
    f=f.figure
    f.savefig('Identified particles in 2nd channel.pdf') 
    frame = pd.Series(np.zeros(len(f1['y'])))
    f1=f1.assign(frame=frame.values)
    green_filter_peaks=len(f1.index)/number_of_cells
    
    f2 = tp.locate(img2,int(particle_size),percentile=int((persentile2)))
    plt.figure(figsize=(50,50))
    plt.title('Identified particles in 3rd channel')
    f=tp.annotate(f2,img2)
    f=f.figure
    f.savefig('Identified particles in 3rd channel.pdf') 
    frame = pd.Series(np.ones(len(f2['y'])))
    f2=f2.assign(frame=frame.values)
    blue_filter_peaks=len(f2.index)/number_of_cells
    
    frames=[f1,f2]
    f = pd.concat(frames)
    tracking_space=int(ts)
    t = tp.link_df(f, int(tracking), memory=2)
    t1=tp.filter_stubs(t,2)
    number_of_part=t1['particle'].nunique()
    
    ##########################################################################################################################################################################################################################
    # Cell by cell analysis
    ##########################################################################################################################################################################################################################        
    
    zero=0
  
    for cells in range(len(np.unique(mask_cyto))):
        if cells!=0:
            selected_masks= mask_cyto==int(cells) # Change this value to the specific cell that you need to select
            removed_mask1 = np.einsum('jk,jk->jk',img1, selected_masks) 
            removed_mask2 = np.einsum('jk,jk->jk',img2, selected_masks) 
            f1 = tp.locate(removed_mask1,int(particle_size),percentile=int((persentile)))
            frame = pd.Series(np.zeros(len(f1['y'])))
            f1=f1.assign(frame=frame.values)
            plt.figure(figsize=(10,10))

            
            f2 = tp.locate(removed_mask2,int(particle_size),percentile=int((persentile2)))
            frame = pd.Series(np.ones(len(f2['y'])))
            f2=f2.assign(frame=frame.values)  
 
                        
            if len(f1)==0 or len(f2)==0:
                zero+=1
            else:        
                frames=[f1,f2]
                f = pd.concat(frames)
                tracking_space=int(ts)
                t = tp.link_df(f, int(tracking), memory=2)
                t1=tp.filter_stubs(t,2)
                if t1['particle'].nunique()==0:
                    zero+=1
    global fp          
    clean=number_of_part/number_of_cells               
    fp=number_of_cells-zero
    if int(fp)==0 or int(fp)<0:
        global part
        part=0
    else:
        part=number_of_part/fp
    if fp<0:
        fp=0
    img=img+1
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size = 25)
    text='IMAGE '+str(img)
    pdf.cell(200, 10, txt = text, ln = 140, align = 'C')
    pdf.output("GFG.pdf") 
    pdf_output(img)
            
    number.append(number_of_cells)
    fp2.append(fp)
    number_of_parti.append(part)
    number_of_part_in_Green.append(green_filter_peaks)
    number_of_part_in_Blue.append(blue_filter_peaks)
    inpu_cell_size.append(sizes[image])
    inpu_part_size.append(particle_sizes[image])
    negative_cells.append(zero)
    clean_div.append(clean)     
    return (number_of_cells,part,green_filter_peaks,blue_filter_peaks,zero,clean,fp)      
   
def tracking2(img1,img2,number_of_cells,mask_cyto,image,img):
    import matplotlib.pyplot as plt
    particle_size=particle_sizes[int(image)]
    ts=particle_sizes[int(image)]
    f1 = tp.locate(img1,int(particle_size),percentile=int((persentile)))
    plt.figure(figsize=(50,50))
    plt.title('Identified particles in 2nd channel')
    f=tp.annotate(f1, img1)
    f=f.figure
    f.savefig('Identified particles in 2nd channel.pdf') 
    frame = pd.Series(np.zeros(len(f1['y'])))
    f1=f1.assign(frame=frame.values)
    green_filter_peaks=len(f1.index)/number_of_cells
    
    f2 = tp.locate(img2,int(particle_size),percentile=int((persentile2)))
    plt.figure(figsize=(50,50))
    plt.title('Identified particles in 3rd channel')
    f=tp.annotate(f2,img2)
    f=f.figure
    f.savefig('Identified particles in 3rd channel.pdf') 
    frame = pd.Series(np.ones(len(f2['y'])))
    f2=f2.assign(frame=frame.values)
    blue_filter_peaks=len(f2.index)/number_of_cells
    
    frames=[f1,f2]
    f = pd.concat(frames)
    tracking_space=int(ts)
    t = tp.link_df(f, int(tracking), memory=2)
    t1=tp.filter_stubs(t,2)
    number_of_part=t1['particle'].nunique()
    
    ##########################################################################################################################################################################################################################
    # Cell by cell analysis
    ##########################################################################################################################################################################################################################        
    
    zero=0
  
    for cells in range(len(np.unique(mask_cyto))):
        if cells!=0:
            selected_masks= mask_cyto==int(cells) # Change this value to the specific cell that you need to select
            removed_mask1 = np.einsum('jk,jk->jk',img1, selected_masks) 
            removed_mask2 = np.einsum('jk,jk->jk',img2, selected_masks) 
            f1 = tp.locate(removed_mask1,int(particle_size),percentile=int((persentile)))
            frame = pd.Series(np.zeros(len(f1['y'])))
            f1=f1.assign(frame=frame.values)
            plt.figure(figsize=(10,10))
            f=tp.annotate(f1, removed_mask1)
            
            f2 = tp.locate(removed_mask2,int(particle_size),percentile=int((persentile2)))
            frame = pd.Series(np.ones(len(f2['y'])))
            f2=f2.assign(frame=frame.values)  
 
                        
            if len(f1)==0 or len(f2)==0:
                zero+=1
            else:        
                frames=[f1,f2]
                f = pd.concat(frames)
                tracking_space=int(ts)
                t = tp.link_df(f, int(tracking), memory=2)
                t1=tp.filter_stubs(t,2)
                if t1['particle'].nunique()==0:
                    zero+=1
    global fp          
    clean=number_of_part/number_of_cells               
    fp=number_of_cells-zero
    if int(fp)==0 or int(fp)<0:
        global part
        part=0
    else:
        part=number_of_part/fp
    if fp<0:
        fp=0
    img=img+1
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size = 25)
    text='IMAGE '+str(img)
    pdf.cell(200, 10, txt = text, ln = 140, align = 'C')
    pdf.output("GFG.pdf") 
    pdf_output2(img)
            
    number.append(number_of_cells)
    fp2.append(fp)
    number_of_parti.append(part)
    number_of_part_in_Green.append(green_filter_peaks)
    number_of_part_in_Blue.append(blue_filter_peaks)
    inpu_cell_size.append(sizes[image])
    inpu_part_size.append(particle_sizes[image])
    negative_cells.append(zero)
    clean_div.append(clean)     
    return (number_of_cells,part,green_filter_peaks,blue_filter_peaks,zero,clean,fp)     

def image_analysis(event=None):
    global i 
    counter=0
    for i in range(len(img_list)):   
        import matplotlib.pyplot as plt
        if i > int(fromimage)-1 and i<int(untilimage)+1:
            img=int(i)-1
            image.append(img+1)
            img_0 = file.get_image(int(img))
            img_0.get_frame(z=0, t=0, c=0)

            frame_list   = [i for i in img_0.get_iter_t(c=0, z=0)] 
            channel_list = [i for i in img_0.get_iter_c(t=0, z=0)] 
            z_list_0     = [i for i in img_0.get_iter_z(t=0, c=int(ch0)-1)]
            z_list_c     = [i for i in img_0.get_iter_z(t=0, c=int(ch1)-1)]
            z_list_2     = [i for i in img_0.get_iter_z(t=0, c=int(ch1)-1)] 
            z_list_3     = [i for i in img_0.get_iter_z(t=0, c=int(ch2)-1)] 
                        
            images=stacking_step(z_list_0,z_list_2,z_list_3)
            mask=masking(saturated_images(z_list_0,z_list_3)[0],saturated_images(z_list_0,z_list_3)[1],counter)
            image1=WBNS(images[1])
            image2=WBNS(images[2])
            removed_mask1 = np.einsum('jk,jk->jk',image1[0], mask[0]) # remove all the background  of the cells 
            removed_mask2 = np.einsum('jk,jk->jk',image2[0], mask[0]) 
            removed_mask_1 = np.einsum('jk,jk->jk',removed_mask1, mask[3]) # remove the nucleus for  noise 
            removed_mask_2 = np.einsum('jk,jk->jk',removed_mask2, mask[3]) 
            img_copy1 = removed_mask_1.copy() # making a copy of our img
            img_gaussian_filter_simga_1 = gaussian_filter(img_copy1, sigma=int((sigm)))
            plotting3(img_gaussian_filter_simga_1,'Greens','Filtered channel 2 image')
            img_copy2 = removed_mask_2.copy() # making a copy of our img
            img_gaussian_filter_simga_2 = gaussian_filter(img_copy2, sigma=int((sigm)))
            import matplotlib.pyplot as plt
            matplotlib.use('agg')
            plotting4(img_gaussian_filter_simga_2,'Blues','Filtered channel 3 image') 
            number_of_cells,part,green_filter_peaks,blue_filter_peaks,zero,clean,fp=tracking3(img_gaussian_filter_simga_1,img_gaussian_filter_simga_2,quantitative_analysis(mask[2])[1],mask[1],counter,img)          
            counter+=1
            
    excel_output(image,number,number_of_parti,number_of_part_in_Green,number_of_part_in_Blue,inpu_cell_size,inpu_part_size, negative_cells, clean_div,fp2)
                        
def image_analysis2(event=None):    
    global i 
    counter=0
    for i in range(len(img_list)):        
        if i > int(fromimage)-1 and i<int(untilimage)+1:             
            img=int(i)-1
            image.append(img+1)
            img_0 = file.get_image(int(img))
            img_0.get_frame(z=0, t=0, c=0)            
            #c1=int(nucleus_channel)-1

            frame_list   = [i for i in img_0.get_iter_t(c=0, z=0)] 
            channel_list = [i for i in img_0.get_iter_c(t=0, z=0)] 
            z_list_0     = [i for i in img_0.get_iter_z(t=0, c=int(ch0)-1)]
            z_list_1     = [i for i in img_0.get_iter_z(t=0, c=int(ch1)-1)]
            z_list_2     = [i for i in img_0.get_iter_z(t=0, c=int(ch2)-1)] 
            images=stacking_step(z_list_0,z_list_1,z_list_2)   
            mask=masking2(saturated_images(z_list_0,z_list_0)[0],saturated_images(z_list_0,z_list_0)[1],counter)
            image1=WBNS(images[1])
            image2=WBNS(images[2])
            removed_mask1 = np.einsum('jk,jk->jk',image1[0], mask[0]) # remove all the background  of the cells 
            removed_mask2 = np.einsum('jk,jk->jk',image2[0], mask[0])      
            img_copy1 = removed_mask1.copy() # making a copy of our img
            img_gaussian_filter_simga_1 = gaussian_filter(img_copy1, sigma=int((sigm)))
            plotting(img_gaussian_filter_simga_1,'Greens','Filtered channel 2 image')
            img_copy2 = removed_mask2.copy() # making a copy of our img
            img_gaussian_filter_simga_2 = gaussian_filter(img_copy2, sigma=int((sigm)))
            plotting(img_gaussian_filter_simga_2,'Blues','Filtered channel 3 image')
            import matplotlib.pyplot as plt
            matplotlib.use('agg')
            number_of_cells,part,green_filter_peaks,blue_filter_peaks,zero,clean,fp=tracking2(img_gaussian_filter_simga_1,img_gaussian_filter_simga_2,quantitative_analysis(mask[2])[1],mask[2],counter,img)          
            counter+=1
    excel_output(image,number,number_of_parti,number_of_part_in_Green,number_of_part_in_Blue,inpu_cell_size,inpu_part_size, negative_cells, clean_div,fp2)
    
def secondary_param(event=None):
    line = Frame(root, height=3, width=406, relief='groove',background='#EA6060')
    line.place(x=50, y=450)
    line = Frame(root, height=3, width=406, relief='groove',background='#EA6060')
    line.place(x=50, y=580)
    line = Frame(root, height=127, width=3, relief='groove',background='#EA6060')
    line.place(x=50, y=453)
    line = Frame(root, height=127, width=3, relief='groove',background='#EA6060')
    line.place(x=453, y=453)
    line = Frame(root, height=16, width=3, relief='groove',bg='#EA6060')
    line.place(x=62, y=474)
    e3 = Entry(root, bd =1, bg='white',textvariable=ps,width=15,highlightbackground='#1A507B',highlightthickness=5)
    e3.place(x=245, y=470)
    intvar3 = IntVar()
    
    def process3(event=None):
        content = e3.get()  
        Ps=ps.get()
        persentil.append(Ps)
        intvar3.set(100)
        line = Frame(root, height=16, width=5, relief='flat',bg='#142841')
        line.place(x=61, y=474)
        l = Label(root, text = content ,bg = "#585E63",fg='#152635',width=4)
        l.config(font =("arial bold", 10))
        l.place(x=392,y=473)
        popupmsg('Suggested sigma level: 2', 'Information')
        
    e3.bind('<Return>', process3)
    e3.wait_variable(intvar3)
    line = Frame(root, height=16, width=3, relief='groove',bg='#EA6060')
    line.place(x=62, y=514)
    
    intvar4 = IntVar()    
    def test1():
        s=var_s.get()
        gaus_f.append(s)
        intvar4.set(100)
        line = Frame(root, height=16, width=3, relief='flat',bg='#142841')
        line.place(x=62, y=514)
        line = Frame(root, height=16, width=3, relief='groove',bg='#EA6060')
        line.place(x=62, y=543)
        popupmsg('Suggested noise level: 2', 'Information')
        
    c1=Checkbutton(root, text="1",variable=var_s,onvalue = 1, offvalue = 0,bg='#142841',command = test1)
    c1.place(x=130, y=510)
    c1=Checkbutton(root, text="2",variable=var_s,onvalue = 2, offvalue = 0,bg='#142841',command = test1)
    c1.place(x=178, y=510)
    c1=Checkbutton(root, text="3",variable=var_s,onvalue = 3, offvalue = 0,bg='#142841',command = test1)
    c1.place(x=223, y=510)
    c1=Checkbutton(root, text="4",variable=var_s,onvalue = 4, offvalue = 0,bg='#142841',command = test1)
    c1.place(x=268, y=510)
    c1=Checkbutton(root, text="5",variable=var_s,onvalue = 5, offvalue = 0,bg='#142841',command = test1)
    c1.place(x=313, y=510)
    c1=Checkbutton(root, text="6",variable=var_s,onvalue = 6, offvalue = 0,bg='#142841',command = test1)
    c1.place(x=359, y=510)
    c1=Checkbutton(root, text="7",variable=var_s,onvalue = 7, offvalue = 0,bg='#142841',command = test1)
    c1.place(x=410, y=510) 
                   
    intvar5 = IntVar()    
    def test2():
        n=var.get() 
        noise_l.append(n)
        line = Frame(root, height=3, width=406, relief='groove',background='#142841')
        line.place(x=50, y=450)
        line = Frame(root, height=3, width=406, relief='groove',background='#142841')
        line.place(x=50, y=580)
        line = Frame(root, height=127, width=3, relief='groove',background='#142841')
        line.place(x=50, y=453)
        line = Frame(root, height=127, width=3, relief='groove',background='#142841')
        line.place(x=453, y=453)
        line = Frame(root, height=16, width=3, relief='groove',bg='#142841')
        line.place(x=62, y=474)
        line = Frame(root, height=16, width=3, relief='groove',bg='#142841')
        line.place(x=62, y=543)
        intvar5.set(100)
        popupmsg1('Ready to start, this might take a while!', 'Initiate')
                   
    d1=Checkbutton(root, text="1",variable=var,onvalue = 1, offvalue = 0,bg='#142841',command = test2)
    d1.place(x=130, y=540)###
    d1=Checkbutton(root, text="2",variable=var,onvalue = 2, offvalue = 0,bg='#142841',command = test2)
    d1.place(x=178, y=540) 
    d1=Checkbutton(root, text="3",variable=var,onvalue = 3, offvalue = 0,bg='#142841',command = test2)
    d1.place(x=223, y=540)
    d1=Checkbutton(root, text="4",variable=var,onvalue = 4, offvalue = 0,bg='#142841',command = test2)
    d1.place(x=268, y=540)
    d1=Checkbutton(root, text="5",variable=var,onvalue = 5, offvalue = 0,bg='#142841',command = test2)
    d1.place(x=313, y=540)
    d1=Checkbutton(root, text="6",variable=var,onvalue = 6, offvalue = 0,bg='#142841',command = test2)
    d1.place(x=359, y=540)
    d1=Checkbutton(root, text="7",variable=var,onvalue = 7, offvalue = 0,bg='#142841',command = test2)
    d1.place(x=410, y=540)
    d1.wait_variable(intvar5)
    
def image_params():
    global img ,img_list,file,first_page,last_page
    #img_path = filedialog.askopenfilename(initialdir=os.getcwd()) 
    file = LifFile(str(img_path))
    img_list = [i for i in file.get_iter_image()] 
    enum=[]
   # fp=IntVar()
    ei = Entry(root, bd =1, bg="#C3CBCE",textvariable=fp,width=4,highlightbackground='#EA6060',highlightthickness=5)
    ei.place(x=312, y=185)
    intvar22 = IntVar()
    line = Frame(root, height=19, width=7, relief='flat',background='#263A55')
    line.place(x=60, y=142)
    def process4(event=None):
        global first_page
        first_page=fp.get() 
        ei = Entry(root, bd =2, bg='#263A55',fg='white',relief='flat',textvariable=fp,width=4,highlightbackground='#263A55',highlightthickness=5,highlightcolor='#263A55')
        ei.place(x=312, y=185)
        eii = Entry(root, bd =1, bg="#C3CBCE",textvariable=lp,width=4,highlightbackground='#EA6060',highlightthickness=5)
        eii.place(x=392, y=185)    
        intvar22.set(100)
    ei.bind('<Return>', process4)
    ei.wait_variable(intvar22)  
  #  lp=IntVar()
    eii = Entry(root, bd =1, bg="#C3CBCE",textvariable=lp,width=4,highlightbackground='#EA6060',highlightthickness=5)
    eii.place(x=392, y=185) 
    intvar33 = IntVar()
    def process5(event=None):
        global last_page
        last_page=lp.get() 
        eii = Entry(root, bd =2, bg='#263A55',fg='white',relief='flat',textvariable=lp,width=4,highlightbackground='#263A55',highlightthickness=5,highlightcolor='#263A55')
        eii.place(x=392, y=185)  
        intvar33.set(100)
    eii.bind('<Return>', process5)
    eii.wait_variable(intvar33)

    for i in range(len(img_list)): 
        if i > int(first_page)-1 and i<int(last_page)+1:
            img=int(i)-1
            enum.append(i)
            img_0 = file.get_image(int(img))
            img_0.get_frame(z=0, t=0, c=0)
            
            c1=int(nucleus_channel)
            c2=int(channel1)
            c3=int(channel2)
            #print(int(c1),print(c2),print(c3))
            frame_list   = [i for i in img_0.get_iter_t(c=0, z=0)] 
            channel_list = [i for i in img_0.get_iter_c(t=0, z=0)] 
            z_list_0     = [i for i in img_0.get_iter_z(t=0, c=int(nucleus_channel)-1)]
            z_list_c     = [i for i in img_0.get_iter_z(t=0, c=int(channel1)-1)]
            z_list_2     = [i for i in img_0.get_iter_z(t=0, c=int(ch_1))] 
            z_list_3     = [i for i in img_0.get_iter_z(t=0, c=int(ch_2))] 
            ns = IntVar()
            ts = IntVar()
            plotting(saturated_images(z_list_0, z_list_c)[0],'Reds','Saturated image: '+ str(img+1))
            plotting3(stacking_step(z_list_0,z_list_2,z_list_3)[1],'Blues','Channel 2 Image: '+ str(img+1))   
            plotting4(stacking_step(z_list_0,z_list_2,z_list_3)[2],'Greens','Channel 3 Image: '+ str(img+1))   
            line = Frame(root, height=16, width=3, relief='groove',bg='#EA6060')
            line.place(x=62, y=344)
            e1 =Entry(root, bd =1, bg='white',textvariable=ns,width=6,highlightbackground='#1A507B',highlightthickness=5)
            e1.place(x=320, y=340)              
            intvar = IntVar()
            def process1(event=None):
                content = e1.get()  
                Ns=ns.get() 
                intvar.set(100)
                l = Label(root, text = content,bg = "#C3CBCE",fg='#152635',width=4)
                l.config(font =("arial bold", 10))
                sizes.append(Ns)
                l.place(x=392,y=343)   
            e1.bind('<Return>', process1)
            e1.wait_variable(intvar)
            line = Frame(root, height=16, width=3, relief='groove',bg='white')
            line.place(x=62, y=344)
            line = Frame(root, height=16, width=3, relief='groove',bg='#EA6060')
            line.place(x=62, y=392)
            e2 = Entry(root, bd =1, bg='white',textvariable=ts,width=6,highlightbackground='#1A507B',highlightthickness=5)
            e2.place(x=320, y=390)
            intvar2 = IntVar()
            def process2(event=None):
                content = e2.get()  
                Ts=ts.get()    
                if Ts % 2:
                    intvar2.set(100)
                    particle_sizes.append(Ts)
                    l = Label(root, text = content,bg = "#C3CBCE",fg='#152635',width=4)
                    l.config(font =("arial bold", 10))
                    l.place(x=392,y=394)
                    import matplotlib.pyplot as plt
                    ploting(saturated_images(z_list_0, z_list_2)[0],'Reds','Saturated image: '+ str(img+1))
                    T = Text(root2, height = 5, width = 52)  
                    l = Label(root2, text = "Diamater input: "+str(ns.get())+" pixels",bg = '#4B637F',fg='#AFB9BF')
                    l.config(font =("arial bold", 8))
                    l.place(x=11,y=202)
                    l = Label(root2, text = "Particle input: "+str(ts.get())+" pixels",bg = '#4B637F',fg='#AFB9BF')
                    l.config(font =("arial bold", 8))
                    l.place(x=11,y=222)
                else:
                    messagebox.showinfo("Task failed succesfully!",  "Please instert only odd numbers (ex:1,3)") 
            e2.bind('<Return>', process2)
            e2.wait_variable(intvar2)
            line = Frame(root, height=16, width=3, relief='groove',bg='white')
            line.place(x=62, y=392)
    secondary_param()

def image_params2():
    global img ,img_list,file,first_page,last_page
    #img_path = filedialog.askopenfilename(initialdir=os.getcwd()) 
    file = LifFile(str(img_path))
    img_list = [i for i in file.get_iter_image()] 
    enum=[]
   # fp=IntVar()
    ei = Entry(root, bd =1, bg="#C3CBCE",textvariable=fp,width=4,highlightbackground='#EA6060',highlightthickness=5)
    ei.place(x=312, y=185)
    intvar22 = IntVar()
    line = Frame(root, height=19, width=7, relief='flat',background='#263A55')
    line.place(x=60, y=142)
    def process4(event=None):
        global first_page
        first_page=fp.get() 
        ei = Entry(root, bd =2, bg='#263A55',fg='white',relief='flat',textvariable=fp,width=4,highlightbackground='#263A55',highlightthickness=5,highlightcolor='#263A55')
        ei.place(x=312, y=185)
        eii = Entry(root, bd =1, bg="#C3CBCE",textvariable=lp,width=4,highlightbackground='#EA6060',highlightthickness=5)
        eii.place(x=392, y=185)    
        intvar22.set(100)
    ei.bind('<Return>', process4)
    ei.wait_variable(intvar22)  
  #  lp=IntVar()
    eii = Entry(root, bd =1, bg="#C3CBCE",textvariable=lp,width=4,highlightbackground='#EA6060',highlightthickness=5)
    eii.place(x=392, y=185) 
    intvar33 = IntVar()
    def process5(event=None):
        global last_page
        last_page=lp.get() 
        eii = Entry(root, bd =2, bg='#263A55',fg='white',relief='flat',textvariable=lp,width=4,highlightbackground='#263A55',highlightthickness=5,highlightcolor='#263A55')
        eii.place(x=392, y=185)  
        intvar33.set(100)
    eii.bind('<Return>', process5)
    eii.wait_variable(intvar33)

    for i in range(len(img_list)): 
        if i > int(first_page)-1 and i<int(last_page)+1:
            img=int(i)-1
            enum.append(i)
            img_0 = file.get_image(int(img))
            img_0.get_frame(z=0, t=0, c=0)
            c1=int(nucleus_channel)
            #print(int(c1),print(c2),print(c3))
            frame_list   = [i for i in img_0.get_iter_t(c=0, z=0)] 
            channel_list = [i for i in img_0.get_iter_c(t=0, z=0)] 
            z_list_0     = [i for i in img_0.get_iter_z(t=0, c=int(nucleus_channel)-1)]
            z_list_2     = [i for i in img_0.get_iter_z(t=0, c=int(ch_1))] 
            z_list_3     = [i for i in img_0.get_iter_z(t=0, c=int(ch_2))] 
            ns = IntVar()
            ts = IntVar()
            plotting(saturated_images(z_list_0, z_list_0)[0],'Reds','Saturated image: '+ str(img+1))
            plotting3(stacking_step(z_list_0,z_list_2,z_list_3)[1],'Blues','Channel 2 Image: '+ str(img+1))   
            plotting4(stacking_step(z_list_0,z_list_2,z_list_3)[2],'Greens','Channel 3 Image: '+ str(img+1))   
            line = Frame(root, height=16, width=3, relief='groove',bg='#EA6060')
            line.place(x=62, y=344)
            e1 =Entry(root, bd =1, bg='white',textvariable=ns,width=6,highlightbackground='#1A507B',highlightthickness=5)
            e1.place(x=320, y=340)              
            intvar = IntVar()
            def process1(event=None):
                content = e1.get()  
                Ns=ns.get() 
                intvar.set(100)
                l = Label(root, text = content,bg = "#C3CBCE",fg='#152635',width=4)
                l.config(font =("arial bold", 10))
                sizes.append(Ns)
                l.place(x=392,y=343)   
            e1.bind('<Return>', process1)
            e1.wait_variable(intvar)
            line = Frame(root, height=16, width=3, relief='groove',bg='white')
            line.place(x=62, y=344)
            line = Frame(root, height=16, width=3, relief='groove',bg='#EA6060')
            line.place(x=62, y=392)
            e2 = Entry(root, bd =1, bg='white',textvariable=ts,width=6,highlightbackground='#1A507B',highlightthickness=5)
            e2.place(x=320, y=390)
            intvar2 = IntVar()
            def process2(event=None):
                content = e2.get()  
                Ts=ts.get()    
                if Ts % 2:
                    intvar2.set(100)
                    particle_sizes.append(Ts)
                    l = Label(root, text = content,bg = "#C3CBCE",fg='#152635',width=4)
                    l.config(font =("arial bold", 10))
                    l.place(x=392,y=394)
                    import matplotlib.pyplot as plt
                    ploting(saturated_images(z_list_0, z_list_2)[0],'Reds','Saturated image: '+ str(img+1))
                    T = Text(root2, height = 5, width = 52)  
                    l = Label(root2, text = "Diamater input: "+str(ns.get())+" pixels",bg = '#4B637F',fg='#AFB9BF')
                    l.config(font =("arial bold", 8))
                    l.place(x=11,y=202)
                    l = Label(root2, text = "Particle input: "+str(ts.get())+" pixels",bg = '#4B637F',fg='#AFB9BF')
                    l.config(font =("arial bold", 8))
                    l.place(x=11,y=222)
                else:
                    messagebox.showinfo("Task failed succesfully!",  "Please instert only odd numbers (ex:1,3)") 
            e2.bind('<Return>', process2)
            e2.wait_variable(intvar2)
            line = Frame(root, height=16, width=3, relief='groove',bg='white')
            line.place(x=62, y=392)
    secondary_param()

def excel_output(image,number,number_of_part,number_of_part_in_Green,number_of_part_in_Blue,inpu_cell_size,inpu_part_size,zero_cell,clean_cell,fp2):
    global name
    t=str(outputfile)+'.xlsx'
    image=pd.DataFrame(image,columns=['Image'])
    number=pd.DataFrame(number,columns=['Number of Totally identified cells'])
    fp2=pd.DataFrame(fp2,columns=['Number of positive cells'])
    number_of_part=pd.DataFrame(number_of_part,columns=['Particles per positive cells with colocalization '])
    number_of_part_in_Green=pd.DataFrame(number_of_part_in_Green,columns=['Particles per cell in the Blue Channel'])
    number_of_part_in_Blue=pd.DataFrame(number_of_part_in_Blue,columns=['Particles per cell in the Green Channel'])
    inpu_cell_size=pd.DataFrame(inpu_cell_size,columns=['Selected nucleus size (px)'])
    inpu_part_size=pd.DataFrame(inpu_part_size,columns=['Selected particle size (px)'])
    zero_cell=pd.DataFrame(negative_cells,columns=['Cells with 0 signal'])
    clean_cell=pd.DataFrame(clean_div,columns=['Particles with colocalization per entirety of cells'])
    result = pd.concat([image,inpu_cell_size,inpu_part_size,number,fp2,zero_cell,number_of_part_in_Green,number_of_part_in_Blue,number_of_part,clean_cell], axis=1, sort=False)  
    writer = pd.ExcelWriter(t, engine='xlsxwriter')
    workbook = writer.book
    worksheet = workbook.add_worksheet('Summary')
    writer.sheets['Summary'] = worksheet
    result.to_excel(writer,sheet_name='Summary', index=False, startrow=1, header=False)
    header_format = workbook.add_format({
        'bg_color': '#5F9EA0',  
        'bold': True,           
        'text_wrap': True,
        'valign': 'top',
        'align': 'center',
        'border': 1})
    # Write the column headers with the defined format.
    for col_num, value in enumerate(result.columns.values):
        worksheet.write(0, col_num, value, header_format)
    writer.save()
    final_merging('2')
    



#---------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                    MAIN SCRIPT
#---------------------------------------------------------------------------------------------------------------------------------------------------

def main(argv):

   global sigm
   global tracking
   global persentile
   global persentile2
   global fromimage
   global untilimage
   global nucls
   global parts
   global ch0
   global ch1
   global ch2
   global model
   global nl
   global outputfile
   
 
   sigm=int(2)
   tracking=int(5)
   persentile=int(90)
   persentile2=int(90)
   ch0=int(1)
   ch1=int(2)
   ch2=int(3)
   nl=int(3)
   outputfile=str('result')
   
   try:
      opts, args = getopt.getopt(argv,"hi:o:m:f:u:p:d:t:s:b:g:w:d",["ifile=","ofile=","model=","from=","until=","pc1=","pc2=","trackingsp=","sigma=","ch1=","ch2=","ch0=","nl="])
   except getopt.GetoptError:
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('test.py -i <inputfile> -o <outputfile> -m <cyto/nuclei> -f <from image> -u <until image>  -p <persentile of channel1 >  -d <percentile of channel2 > -t <trackingspace> -s <signa> -w <ch0> -b <ch1> -g <ch2> -d <noise level> / ex: python test.py -i 020721_ch2DDX6.lif -o ll -m cyto -f 2 -u 3 --pc1 90 --pc2 80 -t 5 -s 2 -w 1 -b 2 -g3' )
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
      elif opt in ("-m", "--model"):
         model = arg
      elif opt in ("-f", "--from"):
         fromimage = arg
      elif opt in ("-u", "--until"):
         untilimage = arg
      elif opt in ("-p", "--pc1"):
         persentile = arg
      elif opt in ("-d", "--pc2"):
         persentile2 = arg
      elif opt in ("-t", "--trackingsp"):          
         tracking = arg
      elif opt in ("-s", "--sigma"):
         sigm = arg
      elif opt in ("-w", "--ch0"):
         ch0 = arg     
      elif opt in ("-b", "--ch1"):
         ch1 = arg
      elif opt in ("-g", "--ch2"):
         ch2 = arg  
      elif opt in ("-d", "--nl"):
         nl = arg  
    
         
   print ('Input file is :', inputfile)
   print ('Selected model :', model)
   print('Analyze images from: ',fromimage,', until: ',untilimage)
   print('Persentile of channel 1: ',persentile)
   print('Persentile of channel 2: ',persentile2)
   print('Sigma: ',sigm)
   print('Nucleus channel is: Channel ', ch0 )
   print('Colocalization channels are: Channels ', ch1 ," and ", ch2)
   print('Trackign space: ', tracking)
   print ('Output file: ', outputfile)
   print('Noise level for WBNS: ', nl)

   global file
   file = LifFile(str(inputfile))
   global img_list
   img_list = [i for i in file.get_iter_image()]
   enum=[]
   len(img_list)
   for i in range(len(img_list)):
       if i > int(fromimage)-1 and i < int(untilimage)+1:
           print('Image ',int(i))
           img=int(i)-1
           enum.append(i)        
           img_0 = file.get_image(int(img))
           img_0.get_frame(z=0, t=0, c=0)
           frame_list   = [i for i in img_0.get_iter_t(c=0, z=0)] 
           channel_list = [i for i in img_0.get_iter_c(t=0, z=0)] 
           z_list_0     = [i for i in img_0.get_iter_z(t=0, c=int(ch0)-1)]
           z_list_1     = [i for i in img_0.get_iter_z(t=0, c=int(ch1)-1)]
           z_list_2     = [i for i in img_0.get_iter_z(t=0, c=int(ch2)-1)]            
           plt.figure(figsize=(10,10))
           plt.imshow(saturated_images(z_list_0, z_list_2)[0],cmap='Reds')
           plt.show(block=True)
    
    
           while True:
               cell_size=input('Define the approximate nucleus size: ')
               if not cell_size.isnumeric():
                   print("Sorry, your response is not valid could you please try again!")
                   continue
               else:
                   sizes.append(cell_size)
                   break      
           while True:
               part_sizes=input('Define the approximate  particle size ( odd numbers only ): ')
               if not part_sizes.isnumeric():
                   print("Sorry, your response is not valid could you please try again!")
                   continue
               else:
                   particle_sizes.append(part_sizes)
                   break
        

   if model == str('cyto') :             
       image_analysis()
   elif model == str('nuclei'):
       image_analysis2()
   else :
       print('Please refine the model of choise!')
    
    
   
if __name__ == "__main__":
   main(sys.argv[1:])
