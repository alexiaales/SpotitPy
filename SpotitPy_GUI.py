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

#---------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                  GLOBAL LISTS
#---------------------------------------------------------------------------------------------------------------------------------------------------

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
model = models.Cellpose(gpu=use_GPU, model_type='cyto') 
pdf_files=[]
noise_l=[]
gaus_f=[]
persentil=[]
persentil2=[]


#---------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                  FORMAT GUI
#---------------------------------------------------------------------------------------------------------------------------------------------------

# contrast border thumbnail 
root = Tk()
root.title("  Channel Colocalization Tool")
root.geometry("1040x1040+0+0")#
root.configure(background='#4B637F')
root.attributes('-alpha',0.98)

my_img = tk.PhotoImage(file = r"C:\Users\Alexia\Desktop\bio2.png") 
root.iconphoto(False, my_img)

root2 = Tk()
root2.title("Viewer")
root2.geometry("200x250+1050+0")#
root2.configure(background='#4B637F')
root2.attributes('-alpha',0.95)

root3 = Toplevel()
root3.title("Channel selection")
root3.geometry("200x150+1050+305")#
root3.configure(background='#4B637F')
root3.attributes('-alpha',0.92)
        
T = Text(root3, height = 5, width = 5)  
l = Label(root3, text = "Channels ",bg = '#4B637F',fg='#AFB9BF')
l.config(font =("arial bold ", 10))
l.place(x=106,y=3) 
  
T =Text(root3, height = 5, width = 5)  
l = Label(root3, text = "Nucleus  ",bg = '#4B637F',fg='#AFB9BF')
l.config(font =("arial bold underline", 10))
l.place(x=5,y=30)
        
var_n = IntVar()
var_c2= IntVar()
var_c3= IntVar()
var_c4= IntVar()


def test1():
    global nucleus_channel
    nucleus_channel=var_n.get()
c1=Checkbutton(root3, text="1",variable=var_n,onvalue = 1,bg='#4B637F',command = test1).place(x=85, y=30)
c2=Checkbutton(root3, text="2",variable=var_n,onvalue = 2,bg='#4B637F',command = test1).place(x=123, y=30)
c3=Checkbutton(root3, text="3",variable=var_n,onvalue = 3,bg='#4B637F',command = test1).place(x=160, y=30)
    
T = Text(root3, height = 5, width = 5)  
l = Label(root3, text = "Cutoplasm  ",bg = '#4B637F',fg='#AFB9BF')
l.config(font =("arial bold underline", 10))
l.place(x=5,y=55)
    
var_c1= IntVar()

def test2():
    global channel1
    channel1=var_c1.get()    
d1=Checkbutton(root3, text="1",variable=var_c1,onvalue = 1,bg='#4B637F',command = test2).place(x=85, y=55)
d2=Checkbutton(root3, text="2",variable=var_c1,onvalue = 2,bg='#4B637F',command = test2).place(x=123, y=55)
d3=Checkbutton(root3, text="3",variable=var_c1,onvalue = 3,bg='#4B637F',command = test2).place(x=160, y=55)
    
T = Text(root3, height = 5, width = 5)  
l = Label(root3, text = "Colocalized  ",bg = '#4B637F',fg='#AFB9BF')
l.config(font =("arial bold underline", 10))
l.place(x=5,y=80)

def test3():
    global channel2
    channel2=var_c2.get()    
    
def test4():
    global channel3 

    channel3=var_c3.get()    

def test5():
    global channel4
    channel4=var_c4.get()    
    
da1=Checkbutton(root3, text="1",variable=var_c2,onvalue = 1,bg='#4B637F',command = test3).place(x=85, y=80)
da2=Checkbutton(root3, text="2",variable=var_c3,onvalue = 2,bg='#4B637F',command = test4).place(x=123, y=80)
da3=Checkbutton(root3, text="3",variable=var_c4,onvalue = 3,bg='#4B637F',command = test5).place(x=160, y=80)   


def on():  
    line = Frame(root3, height=3, width=454, relief='groove',background='#4B637F')
    line.place(x=0, y=1)
    line = Frame(root3, height=3, width=454, relief='groove',background='#4B637F')
    line.place(x=0, y=146)
    line = Frame(root3, height=165, width=3, relief='groove',background='#4B637F')
    line.place(x=1, y=1)
    line = Frame(root3, height=165, width=3, relief='groove',background='#4B637F')
    line.place(x=196, y=1)
btn33 = Button(root3, text="Done", width=20, bg='#142841', fg='white', font=('ariel 11 bold'), relief=GROOVE, command=on)
btn33.place(x=5, y=110)    


def impo(event=None):
    import matplotlib
    import matplotlib.pyplot as plt    
    
def plotting_tocheck():
    matplotlib.use('TkAgg') # <-- THIS MAKES IT FAST!
    plt.figure(figsize=(10,10))
    plt.imshow(z_0,cmap='Reds')
    plt.show()
    impo()
    
def plotting_tocheck2():
    matplotlib.use('TkAgg') # <-- THIS MAKES IT FAST!
    plt.figure(figsize=(10,10))
    plt.imshow(z_1,cmap='Blues')
    plt.show()   
    impo()
    
def plotting_tocheck3():
    matplotlib.use('TkAgg') # <-- THIS MAKES IT FAST!
    plt.figure(figsize=(10,10))
    plt.imshow(z_2,cmap='Greens')
    plt.show()
    impo()
        
root4 = Tk()
root4.title("Channel selection")
root4.geometry("200x138+1050+510")
root4.configure(background='#4B637F')
root4.attributes('-alpha',0.92)

line = Frame(root4, height=3, width=20, relief='groove',background='#142841')
line.place(x=14, y=25) 
btn1 = Button(root4, text="View Channel1", width=18,height=1,bg='#57CAC8', fg='WHITE',font=('ariel 9 bold'), relief=GROOVE, command=plotting_tocheck)
btn1.place(x=50, y=12)
line = Frame(root4, height=3, width=20, relief='groove',background='#142841')
line.place(x=14, y=62)
line = Frame(root4, height=3, width=20, relief='groove',background='#142841')
line.place(x=14, y=68)        
btn1 = Button(root4, text="View Channel2", width=18,height=1,bg='#57CAC8', fg='WHITE',font=('ariel 9 bold'), relief=GROOVE, command=plotting_tocheck2)
btn1.place(x=50, y=55)
line = Frame(root4, height=3, width=20, relief='groove',background='#142841')
line.place(x=14, y=103)
line = Frame(root4, height=3, width=20, relief='groove',background='#142841')
line.place(x=14, y=109)
line = Frame(root4, height=3, width=20, relief='groove',background='#142841')
line.place(x=14, y=115)        
btn1 = Button(root4, text="View Channel3", width=18,height=1,bg='#57CAC8', fg='WHITE',font=('ariel 9 bold'), relief=GROOVE, command=plotting_tocheck3)
btn1.place(x=50, y=98)

   
#---------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                  CANVASES
#---------------------------------------------------------------------------------------------------------------------------------------------------

canvas22 = Canvas(root2, width="165", height= "150", relief=RIDGE, bd=1, bg='white',highlightbackground='#1A507B',highlightthickness=5)
canvas22.place(x=10, y=28)
canvas2 = Canvas(root, width="447", height= "490", relief=RIDGE, bd=1, bg='white',highlightbackground='#1A507B',highlightthickness=5)
canvas2.place(x=540, y=80)
canvas3 = Canvas(root, width="395", height= "490" , bg='#263A55',highlightbackground='#263A55',highlightthickness=5)
canvas3.place(x=50, y=80)
canvas4 = Canvas(root, width="395", height= "120" , bg='#142841',highlightbackground='#142841',highlightthickness=5)
canvas4.place(x=50, y=450)
canvas5 = Canvas(root, width="395", height= "130" , bg='#263A55',highlightbackground='#263A55',highlightthickness=5)
canvas5.place(x=50, y=310)

#---------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                   LINES
#---------------------------------------------------------------------------------------------------------------------------------------------------


line = Frame(root, height=140, width=1, relief='groove')
line.place(x=302, y=310)
line = Frame(root, height=3, width=454, relief='groove',background='#1A507B')
line.place(x=540, y=414)
line = Frame(root, height=165, width=4, relief='groove',background='#1A507B')
line.place(x=770, y=414)
line = Frame(root, height=16, width=3, relief='groove')
line.place(x=62, y=344)
line =Frame(root, height=16, width=3, relief='groove')
line.place(x=62, y=392)
line = Frame(root, height=633, width=1, relief='groove')
line.place(x=490, y=7)
line = Frame(root, height=1, width=405, relief='groove')
line.place(x=50, y=450)
line = Frame(root, height=1, width=405, relief='groove')
line.place(x=50, y=250)
line = Frame(root, height=1, width=405, relief='groove')
line.place(x=50, y=50)
line = Frame(root3, height=1, width=100, relief='groove',background='#142841')
line.place(x=88, y=25)
line = Frame(root, height=140, width=4, relief='groove',background='#142841')
line.place(x=50, y=310)
line = Frame(root, height=140, width=4, relief='groove',background='#142841')
line.place(x=452, y=310)
line = Frame(root, height=137, width=4, relief='groove',background='#142841')
line.place(x=50, y=110)
line = Frame(root, height=137, width=4, relief='groove',background='#142841')
line.place(x=450, y=110)
line = Frame(root, height=4, width=400, relief='groove',background='#142841')
line.place(x=52, y=242)
line = Frame(root2, height=400, width=4, relief='groove',background='#142841')
line.place(x=0, y=0)
line = Frame(root3, height=400, width=4, relief='groove',background='#142841')
line.place(x=0, y=0)
line = Frame(root4, height=400, width=4, relief='groove',background='#142841')
line.place(x=0, y=0)


#---------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                   TEXT
#---------------------------------------------------------------------------------------------------------------------------------------------------

T = Text(root, height = 5, width = 52)  
l = Label(root, text = "IMAGE PREVIEW",bg = "WHITE",fg='#666A68',width=45,height=2)
l.config(font =("helvetica", 12,'bold'))
l.place(x=540,y=20)

l = Label(root, text = "SELECT YOUR FILE",bg = "WHITE",fg='#666A68',width=40,height=2)
l.config(font =("helvetica ", 12, 'bold'))
l.place(x=50,y=20)

l = Label(root, text = "   ",bg = "#57CAC8",fg='#152635',width=4,height=2)
l.config(font =("helvetica ", 12, 'bold'))
l.place(x=50,y=20)

line =Frame(l, height=3, width=23, relief='flat')
line.place(x=9, y=11)
line = Frame(l, height=3, width=23, relief='flat')
line.place(x=9, y=19)
line = Frame(l, height=3, width=23, relief='flat')
line.place(x=9, y=27)

l = Label(root, text = "    ",bg = '#4B637F',fg='WHITE',width=40,height=1)
l.config(font =("helvetica ", 12, 'bold'))
l.place(x=50,y=246)

l = Label(root, text = "SETTING PARAMETERS",bg = "WHITE",fg='#666A68',width=40,height=2)
l.config(font =("helvetica ", 12, 'bold'))
l.place(x=50,y=271)

T = Text(root, height = 5, width = 52)  
l = Label(root, text = "What model would you like: ",bg = '#263A55',fg='#AFB9BF')
l.config(font =("arial bold", 10))
l.place(x=70,y=140)
btt1 = Button(root, text="Cyto", width=6, bg='#57CAC8', fg='white', font=('ariel 9 bold'), relief=GROOVE)
btt1.place(x=310, y=140)
btt2 = Button(root, text="Nuclei", width=6, bg='#57CAC8', fg='white', font=('ariel 9 bold'), relief=GROOVE, command=None)
btt2.place(x=380, y=140)

fp=IntVar()
lp=IntVar()
T = Text(root, height = 5, width = 52)  
l = Label(root, text = "Select images  From: ",bg = '#263A55',fg='#AFB9BF')
l.config(font =("arial bold", 10))
l.place(x=70,y=190)
T = Text(root, height = 5, width = 52)  
l = Label(root, text = "To ",bg = '#263A55',fg='#AFB9BF')
l.config(font =("arial bold", 10))
l.place(x=360,y=190)
ei = Entry(root, bd =1, bg="#C3CBCE",textvariable=fp,width=4,highlightbackground='#1A507B',highlightthickness=5)
ei.place(x=312, y=185)
eii = Entry(root, bd =1, bg="#C3CBCE",textvariable=lp,width=4,highlightbackground='#1A507B',highlightthickness=5)
eii.place(x=392, y=185)

ns = IntVar()
T = Text(root, height = 5, width = 52)  
l = Label(root, text = "Approxiamte nucleus diameter: ",bg = '#263A55',fg='#AFB9BF')
l.config(font =("arial bold", 10))
l.place(x=70,y=340)
l = Label(root, text = "       ",bg = "#C3CBCE",fg='#152635',width=4)
l.config(font =("arial bold", 10))
l.place(x=392,y=343)
e1 =Entry(root, bd =1, bg='white',textvariable=ns,width=6,highlightbackground='#1A507B',highlightthickness=5)
e1.place(x=320, y=340)              

ts = IntVar()
T = Text(root, height = 5, width = 52)  
l = Label(root, text = "Define the particle diameter: ",bg = '#263A55',fg='#AFB9BF')
l.config(font =("arial bold", 10))
l.place(x=70,y=390)
l = Label(root, text = "       ",bg = "#C3CBCE",fg='#152635',width=4)
l.config(font =("arial bold", 10))
l.place(x=392,y=394)
e2 = Entry(root, bd =1, bg='white',textvariable=ts,width=6,highlightbackground='#1A507B',highlightthickness=5)
e2.place(x=320, y=390)

ps = IntVar()
T = Text(root, height = 5, width = 52)  
l = Label(root, text = "Define 1st persentile: ",bg = '#142841',fg='#AFB9BF')
l.config(font =("arial bold", 10))
l.place(x=70,y=465)
l = Label(root, text = "       ",bg = "#C3CBCE",fg='#152635',width=4)
l.config(font =("arial bold", 10))
l.place(x=392,y=468)
e3 = Entry(root, bd =1, bg='white',textvariable=ps,width=15,highlightbackground='#1A507B',highlightthickness=5)
e3.place(x=245, y=465)

ps2 = IntVar()
T = Text(root, height = 5, width = 52)  
l = Label(root, text = "Define 2nd persentile: ",bg = '#142841',fg='#AFB9BF')
l.config(font =("arial bold", 10))
l.place(x=70,y=505)
l = Label(root, text = "       ",bg = "#C3CBCE",fg='#152635',width=4)
l.config(font =("arial bold", 10))
l.place(x=392,y=508)
e4 = Entry(root, bd =1, bg='white',textvariable=ps2,width=15,highlightbackground='#1A507B',highlightthickness=5)
e4.place(x=245, y=505)


def process3(event=None):
    content = e3.get()  
    l = Label(root, text = content ,bg = "#585E63",fg='#152635',width=4)
    global Ps
    Ps=ps.get()
    print(ps.get())
    l.config(font =("arial bold", 10))
    l.place(x=392,y=340)
e3.bind('<Return>', process3)

def process4(event=None):
    content = e4.get()  
    l = Label(root, text = content ,bg = "#585E63",fg='#152635',width=4)
    global Ps
    Ps2=ps2.get()
    print(ps.get())
    l.config(font =("arial bold", 10))
    l.place(x=392,y=440)
e4.bind('<Return>', process4)

T = Text(root2, height = 5, width = 52)  
l = Label(root2, text = "Previous Image: ",bg = '#4B637F',fg='#AFB9BF')
l.config(font =("arial bold", 10))
l.place(x=11,y=5)



#-----------------------------------------------------------------------------------------------------------------------------------------------
#
#                                                                  FUNCTIONS

#-----------------------------------------------------------------------------------------------------------------------------------------------


def model_selection(event=None):
    global img_path,var_model
    pdf_files.append('Intro.pdf') 
    line = Frame(root, height=16, width=3, relief='groove')
    line.place(x=62, y=144)
    img_path = filedialog.askopenfilename(initialdir=os.getcwd())
    line = Frame(root, height=16, width=3, relief='groove',background='#EA6060')
    line.place(x=62, y=144)
    var_model = tk.IntVar()
    btt1 = Button(root, text="Cyto", width=6, bg='#57CAC8', fg='white', font=('ariel 9 bold'), relief=GROOVE, command=lambda: var_model.set(1))
    btt1.place(x=310, y=140)
    btt2 = Button(root, text="Nuclei", width=6, bg='#57CAC8', fg='white', font=('ariel 9 bold'), relief=GROOVE, command=lambda: var_model.set(2))
    btt2.place(x=380, y=140)
    root.wait_variable(var_model) 
    
    if var_model.get()==1:
        btt1 = Button(root, text="Cyto", width=6, bg='#263A55', fg='white', font=('ariel 9 bold'), relief=GROOVE, command=lambda: var_model.set(1))
        btt1.place(x=310, y=140)

    else :    
        btt2 = Button(root, text="Nuclei", width=6, bg='#263A55', fg='white', font=('ariel 9 bold'), relief=GROOVE, command=lambda: var_model.set(2))
        btt2.place(x=380, y=140)
        
    line = Frame(root, height=19, width=7, relief='flat',background='#263A55')
    line.place(x=60, y=142)
    
    if var_model.get()==1:
    
    
        root3 = Toplevel()
        root3.title("Channel selection")
        root3.geometry("200x150+1050+305")#
        root3.configure(background='#4B637F')
        root3.attributes('-alpha',0.99)
            
        T = Text(root3, height = 5, width = 5)  
        l = Label(root3, text = "Channels ",bg = '#4B637F',fg='#AFB9BF')
        l.config(font =("arial bold ", 10))
        l.place(x=106,y=3) 
        
        T = Text(root3, height = 5, width = 5)  
        l = Label(root3, text = "Nucleus  ",bg = '#4B637F',fg='#AFB9BF')
        l.config(font =("arial bold underline", 10))
        l.place(x=5,y=30)
            
        line = Frame(root3, height=1, width=100, relief='groove',background='#142841')
        line.place(x=88, y=25)    
        line = Frame(root3, height=3, width=454, relief='groove',background='#EA6060')
        line.place(x=0, y=1)
        line = Frame(root3, height=3, width=454, relief='groove',background='#EA6060')
        line.place(x=0, y=146)
        line = Frame(root3, height=165, width=3, relief='groove',background='#EA6060')
        line.place(x=1, y=1)
        line = Frame(root3, height=165, width=3, relief='groove',background='#EA6060')
        line.place(x=196, y=1)
        
        var_n = IntVar()
        def test1():
            global nucleus_channel
            nucleus_channel=var_n.get()
        c1=Checkbutton(root3, text="1",variable=var_n,onvalue = 1,bg='#4B637F',command = test1).place(x=85, y=30)
        c2=Checkbutton(root3, text="2",variable=var_n,onvalue = 2,bg='#4B637F',command = test1).place(x=123, y=30)
        c3=Checkbutton(root3, text="3",variable=var_n,onvalue = 3,bg='#4B637F',command = test1).place(x=160, y=30)
            
        T = Text(root3, height = 5, width = 5)  
        l = Label(root3, text = "Cutoplasm  ",bg = '#4B637F',fg='#AFB9BF')
        l.config(font =("arial bold underline", 10))
        l.place(x=5,y=55)
            
        var_c1= IntVar()

  
        def test2():
            global channel1
            channel1=var_c1.get()    
        d1=Checkbutton(root3, text="1",variable=var_c1,onvalue = 1,bg='#4B637F',command = test2).place(x=85, y=55)
        d2=Checkbutton(root3, text="2",variable=var_c1,onvalue = 2,bg='#4B637F',command = test2).place(x=123, y=55)
        d3=Checkbutton(root3, text="3",variable=var_c1,onvalue = 3,bg='#4B637F',command = test2).place(x=160, y=55)
            
        T = Text(root3, height = 5, width = 5)  
        l = Label(root3, text = "Colocalized  ",bg = '#4B637F',fg='#AFB9BF')
        l.config(font =("arial bold underline", 10))
        l.place(x=5,y=80)
        
        
        var_c2= IntVar()
        var_c3= IntVar()
        var_c4= IntVar()   
        
        global channel2
        global channel3
        global channel4
        
        channel2=0
        channel3=0
        channel4=0
        
        def test3():
            global channel2
            channel2=0
            channel2=var_c2.get()    
            
        def test4():
            global channel3 
            channel3=0
            channel3=var_c3.get()    
        
        def test5():
            global channel4
            channel4=0
            channel4=var_c4.get()    
            
        da1=Checkbutton(root3, text="1",variable=var_c2,onvalue = 1,bg='#4B637F',command = test3).place(x=85, y=80)
        da2=Checkbutton(root3, text="2",variable=var_c3,onvalue = 2,bg='#4B637F',command = test4).place(x=123, y=80)
        da3=Checkbutton(root3, text="3",variable=var_c4,onvalue = 3,bg='#4B637F',command = test5).place(x=160, y=80)   
    
        cha=[]
        def on():

            channel4a=int(channel4)-1
            channel3a=int(channel3)-1
            channel2a=int(channel2)-1
            ch=[channel4a,channel3a,channel2a]
      
            for x in ch :
                if x>=0:
                    cha.append(x)
            global ch_1
            global ch_2 

            ch_1=cha[0]
            ch_2=cha[1]
            line = Frame(root3, height=3, width=454, relief='groove',background='#4B637F')
            line.place(x=0, y=1)
            line = Frame(root3, height=3, width=454, relief='groove',background='#4B637F')
            line.place(x=0, y=146)
            line = Frame(root3, height=165, width=3, relief='groove',background='#4B637F')
            line.place(x=1, y=1)
            line = Frame(root3, height=165, width=3, relief='groove',background='#4B637F')
            line.place(x=196, y=1)
            image_params()
            
        btn33 = Button(root3, text="Done", width=20, bg='#142841', fg='white', font=('ariel 11 bold'), relief=GROOVE, command=on)
        btn33.place(x=5, y=106)    
        root3.mainloop()  
        
    else :
        
        root3 = Toplevel()
        root3.title("Channel selectioooon")
        root3.geometry("200x150+1050+305")#
        root3.configure(background='#4B637F')
        root3.attributes('-alpha',0.99)
            
        T = Text(root3, height = 5, width = 5)  
        l = Label(root3, text = "Channels ",bg = '#4B637F',fg='#AFB9BF')
        l.config(font =("arial bold ", 10))
        l.place(x=106,y=5) 
        
        T = Text(root3, height = 5, width = 5)  
        l = Label(root3, text = "Nucleus  ",bg = '#4B637F',fg='#AFB9BF')
        l.config(font =("arial bold underline", 10))
        l.place(x=5,y=40)
        
        
        var_c2= IntVar()
        var_c3= IntVar()
        var_c4= IntVar()

            
        line = Frame(root3, height=1, width=100, relief='groove',background='#142841')
        line.place(x=88, y=30)    
        line = Frame(root3, height=3, width=454, relief='groove',background='#EA6060')
        line.place(x=0, y=1)
        line = Frame(root3, height=3, width=454, relief='groove',background='#EA6060')
        line.place(x=0, y=146)
        line = Frame(root3, height=165, width=3, relief='groove',background='#EA6060')
        line.place(x=1, y=1)
        line = Frame(root3, height=165, width=3, relief='groove',background='#EA6060')
        line.place(x=196, y=1)
        
        var_n = IntVar()
        def test1():
            global nucleus_channel
            nucleus_channel=var_n.get()
        c1=Checkbutton(root3, text="1",variable=var_n,onvalue = 1,bg='#4B637F',command = test1).place(x=85, y=40)
        c2=Checkbutton(root3, text="2",variable=var_n,onvalue = 2,bg='#4B637F',command = test1).place(x=123, y=40)
        c3=Checkbutton(root3, text="3",variable=var_n,onvalue = 3,bg='#4B637F',command = test1).place(x=160, y=40)
                                      
        T = Text(root3, height = 5, width = 5)  
        l = Label(root3, text = "Colocalized  ",bg = '#4B637F',fg='#AFB9BF')
        l.config(font =("arial bold underline", 10))
        l.place(x=5,y=70)
    
               
        global channel22
        global channel32
        global channel42
        channel22=0
        channel32=0
        channel42=0
                
        def test3():
            global channel22
            channel22=0
            channel22=var_c2.get()    
                   
        def test4():
            global channel32 
            channel32=0
            channel32=var_c3.get()    
                
        def test5():
            global channel42
            channel42=0
            channel42=var_c4.get()     
            
        da1=Checkbutton(root3, text="1",variable=var_c2,onvalue = 1,bg='#4B637F',command = test3).place(x=85, y=70)
        da2=Checkbutton(root3, text="2",variable=var_c3,onvalue = 2,bg='#4B637F',command = test4).place(x=123, y=70)
        da3=Checkbutton(root3, text="3",variable=var_c4,onvalue = 3,bg='#4B637F',command = test5).place(x=160, y=70)   
            
        cha=[]
        def on():
            channel4a=int(channel42)-1
            channel3a=int(channel32)-1
            channel2a=int(channel22)-1
            ch=[channel4a,channel3a,channel2a]
    
            for x in ch :
                if x>=0:
                    cha.append(x)
            global ch_1
            global ch_2        
            ch_1=cha[0]
            ch_2=cha[1]
            line = Frame(root3, height=3, width=454, relief='groove',background='#4B637F')
            line.place(x=0, y=1)
            line = Frame(root3, height=3, width=454, relief='groove',background='#4B637F')
            line.place(x=0, y=146)
            line = Frame(root3, height=165, width=3, relief='groove',background='#4B637F')
            line.place(x=1, y=1)
            line = Frame(root3, height=165, width=3, relief='groove',background='#4B637F')
            line.place(x=196, y=1)
            image_params2()
                        
        btn33 = Button(root3, text="Done", width=20, bg='#142841', fg='white', font=('ariel 11 bold'), relief=GROOVE, command=on)
        btn33.place(x=5, y=106)    
        root3.mainloop()  
        
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
    per='Persentile: '+str(persentil)+','+str(persentil2)
    filt='Sigma Gaussian filter: '+str(gaus_f[0])
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
    os.remove(os.path.join('./', 'Saturated image'))
    os.remove(os.path.join('./', 'Channel 2 Image'))
    os.remove(os.path.join('./', 'Channel 3 Image'))
    for pdf in pdf_files:
        merger.append(pdf)
    merger.write(str('Report')+'.pdf')
    merger.close()
    for x in pdf_files:
        os.remove(os.path.join('./', x))
    popupmsg('Analysis completed','Terminal')

def get_going():  
    if var_model.get()==1:
        image_analysis()
    else:
        image_analysis2()

def popupmsg1(msg, title):
    root = tk.Tk()
    root.title(title)
    root.geometry("220x80+550+300")
    
    root.attributes('-alpha',0.80)
    root.configure(background='#581845')
    label = Label(root, text=msg,background='#581845',fg='white')
    label.pack(side="top", fill="x", pady=10)
    B1 = tk.Button(root, text="START",activebackground="#142841", command = lambda:[get_going(),root.destroy()])
    B1.pack()
    root.mainloop()
    
def popupmsg(msg, title):

    root = tk.Tk()
    root.title(title)
    root.geometry("140x80+550+300")#
    root.attributes('-alpha',0.80)
    root.configure(background='#581845')
    label = Label(root, text=msg,background='#581845',fg='white')
    label.pack(side="top", fill="x", pady=10)
    B1 = tk.Button(root, text="Okay", command = root.destroy)
    B1.pack()
      
def plot(image,color,plot_title):     
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig=plt.figure(figsize=(0.8,0.8))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)
    plt.show()       
    plt.figure(figsize=(10,10))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)
    plt.savefig(plot_title+'.pdf') 
    
def ploting(image,color,plot_title):     
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    canvas22 = Canvas(root2, width="165", height= "150", relief=RIDGE, bd=1, bg='white',highlightbackground='#1A507B',highlightthickness=5)
    canvas22.place(x=15, y=37) 
    fig=plt.figure(figsize=(0.8,0.8))
    plt.imshow(image,cmap=color)
    plt.show()       
    canvas = FigureCanvasTkAgg(fig, master=canvas22)
    canvas.get_tk_widget().grid(row=0, column=0, ipadx=55, ipady=40)
    plt.figure(figsize=(10,10))
    plt.imshow(image,cmap=color)
    plt.savefig(plot_title+'.pdf') 
    
def plotting(image,color,plot_title): 
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    canvas2 = Canvas(root, width="450", height= "490", relief=RIDGE, bd=1, bg='white',highlightbackground='#1A507B',highlightthickness=5)
    canvas2.place(x=550, y=86)
    fig=plt.figure(figsize=(5,4))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)
    plt.show()       
    canvas = FigureCanvasTkAgg(fig, master=canvas2)
    canvas.get_tk_widget().grid(row=0, column=0, ipadx=40, ipady=20)
    plt.figure(figsize=(10,10))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)
    plt.savefig(plot_title+'.pdf') 
                
def plotting2(image,color,plot_title): 
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    canvas2 = Canvas(root, width="450", height= "490", relief=RIDGE, bd=1, bg='white',highlightbackground='#1A507B',highlightthickness=5)
    canvas2.place(x=549, y=425)
    fig=plt.figure(figsize=(1.5,1.5))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)
    plt.show()       
    canvas = FigureCanvasTkAgg(fig, master=canvas2)
    canvas.get_tk_widget().grid(row=0, column=0, ipadx=30, ipady=20)  
    plt.figure(figsize=(10,10))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)
    plt.savefig(plot_title+'.pdf')       

def plotting3(image,color,plot_title):
    matplotlib.use('agg')
    import matplotlib.pyplot as plt    
    canvas2 = Canvas(root, width="450", height= "490", relief=RIDGE, bd=1, bg='white',highlightbackground='#1A507B',highlightthickness=5)
    canvas2.place(x=575, y=425)
    fig=plt.figure(figsize=(1.5,1.5))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)    
    plt.show()       
    canvas = FigureCanvasTkAgg(fig, master=canvas2)
    canvas.get_tk_widget().grid(row=0, column=0, ipadx=30, ipady=20)  
    plt.figure(figsize=(10,10))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)
    plt.savefig(plot_title+'.pdf') 

def plotting4(image,color,plot_title):     
    matplotlib.use('agg')
    import matplotlib.pyplot as plt   
    canvas2 = Canvas(root, width="450", height= "490", relief=RIDGE, bd=1, bg='white',highlightbackground='#1A507B',highlightthickness=5)
    canvas2.place(x=801, y=425)
    fig=plt.figure(figsize=(1.5,1.5))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)
    plt.show()       
    canvas = FigureCanvasTkAgg(fig, master=canvas2)
    canvas.get_tk_widget().grid(row=0, column=0, ipadx=30, ipady=20)  
    plt.figure(figsize=(10,10))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)
    plt.savefig(plot_title+'.pdf') 
    
def plotting5(image,color,plot_title):     
    matplotlib.use('agg')
    import matplotlib.pyplot as plt   
    canvas2 = Canvas(root, width="450", height= "490", relief=RIDGE, bd=1, bg='white',highlightbackground='#1A507B',highlightthickness=5)
    canvas2.place(x=801, y=425)
    fig=plt.figure(figsize=(1.5,1.5))
    plt.imshow(image,cmap=color)
    plt.title(plot_title)
    plt.show()       
    canvas = FigureCanvasTkAgg(fig, master=canvas2)
    canvas.get_tk_widget().grid(row=0, column=0, ipadx=30, ipady=20)  
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
    noise_lvl = int(3) #default = 1
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
  # model = models.Cellpose(gpu=use_GPU, model_type='nuclei')  
    masks, flows, styles, diams = model.eval(nucleus, diameter=int(cell_s), flow_threshold=None, channels=[0,0])
    masks= gaussian_filter(masks, sigma=int(gaus_f[0]))
    plotting3(masks,'Reds','Masked nucleus')

  # model = models.Cellpose(gpu=use_GPU, model_type='cyto') 
    masksc, flowsc, stylesc, diamsc = model.eval(cytoplasm, diameter=int(cell_s)*2, flow_threshold=None, channels=[0,0])
    masksc= gaussian_filter(masksc, sigma=int(gaus_f[0]))
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
   # model = models.Cellpose(gpu=use_GPU, model_type='nuclei') 
    cell_s=sizes[int(image)] 
    masks, flows, styles, diams = model.eval(nucleus, diameter=int(cell_s), flow_threshold=None, channels=[0,0])
    masks_copy = masks.copy()
    masks= gaussian_filter(masks, sigma=int(gaus_f[0]))
    plotting(masks,'Oranges','Masked nucleus')
    
    masksc, flowsc, stylesc, diamsc = model.eval(cytoplasm, diameter=int(cell_s)*1.5, flow_threshold=None, channels=[0,0])
    masksc= gaussian_filter(masksc, sigma=int(gaus_f[0]))
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

def tracking(img1,img2,number_of_cells,mask_cyto,image,img):
    import matplotlib.pyplot as plt
    particle_size=particle_sizes[int(image)]
    ts=particle_sizes[int(image)]
    f1 = tp.locate(img1,int(particle_size),percentile=int((persentil)[0]))
    plt.figure(figsize=(50,50))
    plt.title('Identified particles in 2nd channel')
    f=tp.annotate(f1, img1)
    f=f.figure
    f.savefig('Identified particles in 2nd channel.pdf') 
    frame = pd.Series(np.zeros(len(f1['y'])))
    f1=f1.assign(frame=frame.values)
    green_filter_peaks=len(f1.index)/number_of_cells
    
    f2 = tp.locate(img2,int(particle_size),percentile=int((persentil2)[0]))
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
    t = tp.link_df(f, int(tracking_space), memory=2)
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
                    
            f1 = tp.locate(removed_mask1,int(particle_size),percentile=int((persentil)[0]))
            frame = pd.Series(np.zeros(len(f1['y'])))
            f1=f1.assign(frame=frame.values)
            f2 = tp.locate(removed_mask2,int(particle_size),percentile=int((persentil2)[0]))
            frame = pd.Series(np.ones(len(f2['y'])))
            f2=f2.assign(frame=frame.values)        

            if len(f1)==0 or len(f2)==0:
                zero+=1
 
            else:        
                frames=[f1,f2]
                f = pd.concat(frames)
                tracking_space=int(ts)
                t = tp.link_df(f, int(tracking_space), memory=2)
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
    
    particle_size=particle_sizes[int(image)]
    ts=particle_sizes[int(image)]
    f1 = tp.locate(img1,int(particle_size),percentile=int((persentil)[0]))
    plt.figure(figsize=(50,50))
    plt.title('Identified particles in 2nd channel')
    f=tp.annotate(f1, img1)
    f=f.figure
    f.savefig('Identified particles in 2nd channel.pdf') 
    frame = pd.Series(np.zeros(len(f1['y'])))
    f1=f1.assign(frame=frame.values)
    green_filter_peaks=len(f1.index)/number_of_cells
    
    f2 = tp.locate(img2,int(particle_size),percentile=int((persentil2)[0]))
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
    t = tp.link_df(f, int(tracking_space), memory=2)
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
            f1 = tp.locate(removed_mask1,int(particle_size),percentile=int((persentil)[0]))
            frame = pd.Series(np.zeros(len(f1['y'])))
            f1=f1.assign(frame=frame.values)
            plt.figure(figsize=(10,10))
            f=tp.annotate(f1, removed_mask1)
            
            f2 = tp.locate(removed_mask2,int(particle_size),percentile=int((persentil2)[0]))
            frame = pd.Series(np.ones(len(f2['y'])))
            f2=f2.assign(frame=frame.values)  
            plt.figure(figsize=(10,10))
            f=tp.annotate(f2, removed_mask2)
                        
            if len(f1)==0 or len(f2)==0:
                zero+=1
            else:        
                frames=[f1,f2]
                f = pd.concat(frames)
                tracking_space=int(ts)
                t = tp.link_df(f, int(tracking_space), memory=2)
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
    print(number_of_cells,part,green_filter_peaks,blue_filter_peaks,zero,clean,fp)
    
    return (number_of_cells,part,green_filter_peaks,blue_filter_peaks,zero,clean,fp)     
 
def image_analysis(event=None):
    global i 
    counter=0
    for i in range(len(img_list)):   
        import matplotlib.pyplot as plt
        if i > int(first_page)-1 and i<int(last_page)+1:
            img=int(i)-1
            image.append(img+1)
            img_0 = file.get_image(int(img))
            img_0.get_frame(z=0, t=0, c=0)
            c1=int(nucleus_channel)
            c2=int(channel1)
            #print(int(c1),print(c2),print(c3))
            frame_list   = [i for i in img_0.get_iter_t(c=0, z=0)] 
            channel_list = [i for i in img_0.get_iter_c(t=0, z=0)] 
            z_list_0     = [i for i in img_0.get_iter_z(t=0, c=int(nucleus_channel)-1)]
            z_list_c     = [i for i in img_0.get_iter_z(t=0, c=int(channel1)-1)]
            z_list_2     = [i for i in img_0.get_iter_z(t=0, c=int(ch_1))] 
            z_list_3     = [i for i in img_0.get_iter_z(t=0, c=int(ch_2))] 
      
            images=stacking_step(z_list_0,z_list_2,z_list_3)
            mask=masking(saturated_images(z_list_0,z_list_c)[0],saturated_images(z_list_0,z_list_c)[1],counter)
            image1=WBNS(images[1])
            image2=WBNS(images[2])
            removed_mask1 = np.einsum('jk,jk->jk',image1[0], mask[0]) # remove all the background  of the cells 
            removed_mask2 = np.einsum('jk,jk->jk',image2[0], mask[0]) 
            removed_mask_1 = np.einsum('jk,jk->jk',removed_mask1, mask[3]) # remove the nucleus for  noise 
            removed_mask_2 = np.einsum('jk,jk->jk',removed_mask2, mask[3]) 
            img_copy1 = removed_mask_1.copy() # making a copy of our img
            img_gaussian_filter_simga_1 = gaussian_filter(img_copy1, sigma=int((gaus_f)[0]))
            import matplotlib.pyplot as plt
            plotting3(img_gaussian_filter_simga_1,'Greens','Filtered channel 2 image')
            img_copy2 = removed_mask_2.copy() # making a copy of our img
            img_gaussian_filter_simga_2 = gaussian_filter(img_copy2, sigma=int((gaus_f)[0]))
            import matplotlib.pyplot as plt
            plotting4(img_gaussian_filter_simga_2,'Blues','Filtered channel 3 image')      
            number_of_cells,part,green_filter_peaks,blue_filter_peaks,zero,clean,fp=tracking(img_gaussian_filter_simga_1,img_gaussian_filter_simga_2,quantitative_analysis(mask[2])[1],mask[1],counter,img)
            counter+=1   
            import matplotlib.pyplot as plt               
    
    excel_output(image,number,number_of_parti,number_of_part_in_Green,number_of_part_in_Blue,inpu_cell_size,inpu_part_size, negative_cells, clean_div,fp2)
                        
def image_analysis2(event=None):    
    global i 
    counter=0
    for i in range(len(img_list)):        
        if i > int(first_page)-1 and i<int(last_page)+1:
            
            img=int(i)-1
            image.append(img+1)
            img_0 = file.get_image(int(img))
            img_0.get_frame(z=0, t=0, c=0)            
            c1=int(nucleus_channel)-1

            frame_list   = [i for i in img_0.get_iter_t(c=0, z=0)] 
            channel_list = [i for i in img_0.get_iter_c(t=0, z=0)] 
            z_list_0     = [i for i in img_0.get_iter_z(t=0, c=int(nucleus_channel)-1)]
            z_list_1     = [i for i in img_0.get_iter_z(t=0, c=int(ch_1))]
            z_list_2     = [i for i in img_0.get_iter_z(t=0, c=int(ch_2))] 
            images=stacking_step(z_list_0,z_list_1,z_list_2)   
            mask=masking2(saturated_images(z_list_0,z_list_0)[0],saturated_images(z_list_0,z_list_0)[1],counter)
            image1=WBNS(images[1])
            image2=WBNS(images[2])
            removed_mask1 = np.einsum('jk,jk->jk',image1[0], mask[0]) # remove all the background  of the cells 
            removed_mask2 = np.einsum('jk,jk->jk',image2[0], mask[0])      
            img_copy1 = removed_mask1.copy() # making a copy of our img
            img_gaussian_filter_simga_1 = gaussian_filter(img_copy1, sigma=int((gaus_f)[0]))
            plotting(img_gaussian_filter_simga_1,'Greens','Filtered channel 2 image')
            img_copy2 = removed_mask2.copy() # making a copy of our img
            img_gaussian_filter_simga_2 = gaussian_filter(img_copy2, sigma=int((gaus_f)[0]))
            plotting(img_gaussian_filter_simga_2,'Blues','Filtered channel 3 image')
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
    line.place(x=62, y=468)
    e3 = Entry(root, bd =1, bg='white',textvariable=ps,width=15,highlightbackground='#1A507B',highlightthickness=5)
    e3.place(x=245, y=465)
    
    intvar3 = IntVar()
    intvar33= IntVar()
    
    def process3(event=None):
        content = e3.get()  
        Ps=ps.get()
        persentil.append(Ps)
        intvar3.set(100)
        line = Frame(root, height=16, width=5, relief='flat',bg='#142841')
        line.place(x=61, y=468)
        l = Label(root, text = content ,bg = "#585E63",fg='#152635',width=4)
        l.config(font =("arial bold", 10))
        l.place(x=392,y=468)
        
    e3.bind('<Return>', process3)
    e3.wait_variable(intvar3)
    line = Frame(root, height=16, width=3, relief='groove',bg='#EA6060')
    line.place(x=62, y=508)
   
    e4 = Entry(root, bd =1, bg='white',textvariable=ps2,width=15,highlightbackground='#1A507B',highlightthickness=5)
    e4.place(x=245, y=505)
    
    def process4(event=None):
        content2 = e4.get()  
        Ps2=ps2.get()
        persentil2.append(Ps2)
        intvar33.set(100)
        line = Frame(root, height=16, width=3, relief='groove',bg='#142841')
        line.place(x=62, y=508)                
        line = Frame(root, height=16, width=5, relief='flat',bg='#142841')
        line.place(x=61, y=474)
        l = Label(root, text = content2 ,bg = "#585E63",fg='#152635',width=4)
        l.config(font =("arial bold", 10))
        l.place(x=392,y=508)
        
        popupmsg('Suggested sigma level: 2', 'Information')
        
    e4.bind('<Return>', process4)
    e4.wait_variable(intvar33)
    line = Frame(root, height=16, width=3, relief='groove',bg='#EA6060')
    line.place(x=62, y=543)
            
    intvar4 = IntVar()    
    def test1():
        s=var_s.get()
        gaus_f.append(s)
        intvar4.set(100)
        line = Frame(root, height=16, width=3, relief='flat',bg='#142841')
        line.place(x=62, y=514)
        line = Frame(root, height=16, width=3, relief='groove',bg='#EA6060')
        line.place(x=62, y=543)
        line = Frame(root, height=16, width=3, relief='groove',bg='#142841')
        line.place(x=62, y=543)
        popupmsg1('Ready to start, this might take a while!', 'Initiate')
               
    c1=Checkbutton(root, text="1",variable=var_s,onvalue = 1, offvalue = 0,bg='#142841',command = test1)
    c1.place(x=130, y=540)
    c1=Checkbutton(root, text="2",variable=var_s,onvalue = 2, offvalue = 0,bg='#142841',command = test1)
    c1.place(x=178, y=540)
    c1=Checkbutton(root, text="3",variable=var_s,onvalue = 3, offvalue = 0,bg='#142841',command = test1)
    c1.place(x=223, y=540)
    c1=Checkbutton(root, text="4",variable=var_s,onvalue = 4, offvalue = 0,bg='#142841',command = test1)
    c1.place(x=268, y=540)
    c1=Checkbutton(root, text="5",variable=var_s,onvalue = 5, offvalue = 0,bg='#142841',command = test1)
    c1.place(x=313, y=540)
    c1=Checkbutton(root, text="6",variable=var_s,onvalue = 6, offvalue = 0,bg='#142841',command = test1)
    c1.place(x=359, y=540)
    c1=Checkbutton(root, text="7",variable=var_s,onvalue = 7, offvalue = 0,bg='#142841',command = test1)
    c1.place(x=410, y=540) 
    
    noise_l.append(int(3))
                       
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
    t='output.xlsx'
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
    
    
            

#-----------------------------------------------------------------------------------------------------------------------------------------------
#                                                                  CHECK BOXES
#-----------------------------------------------------------------------------------------------------------------------------------------------

T = Text(root, height = 5, width = 52)  
l = Label(root, text = "Sigma : ",bg = '#142841',fg='#AFB9BF')
l.config(font =("arial bold", 10))
l.place(x=70,y=541)
var_s = IntVar()
def test1():
    s=var_s.get()
    var_s = IntVar()
c1=Checkbutton(root, text="1",variable=var_s,onvalue = 1, offvalue = 0,bg='#142841',command = test1).place(x=130, y=540)
c2=Checkbutton(root, text="2",variable=var_s,onvalue = 2, offvalue = 0,bg='#142841',command = test1).place(x=178, y=540)
c3=Checkbutton(root, text="3",variable=var_s,onvalue = 3, offvalue = 0,bg='#142841',command = test1).place(x=223, y=540)
c4=Checkbutton(root, text="4",variable=var_s,onvalue = 4, offvalue = 0,bg='#142841',command = test1).place(x=268, y=540)
c5=Checkbutton(root, text="5",variable=var_s,onvalue = 5, offvalue = 0,bg='#142841',command = test1).place(x=313, y=540)
c6=Checkbutton(root, text="6",variable=var_s,onvalue = 6, offvalue = 0,bg='#142841',command = test1).place(x=359, y=540)
c7=Checkbutton(root, text="7",variable=var_s,onvalue = 7, offvalue = 0,bg='#142841',command = test1).place(x=410, y=540)
#T = Text(root, height = 5, width = 52)  
#l = Label(root, text = "Noise: ",bg = '#142841',fg='#AFB9BF')
#l.config(font =("arial bold", 10))
#l.place(x=70,y=541)
#var= IntVar()
#def test2():
#    n=var.get()    
#var = IntVar()
#d1=Checkbutton(root, text="1",variable=var,onvalue = 1, offvalue = 0,bg='#142841',command = test2).place(x=130, y=540)
#d2=Checkbutton(root, text="2",variable=var,onvalue = 2, offvalue = 0,bg='#142841',command = test2).place(x=178, y=540)
#d3=Checkbutton(root, text="3",variable=var,onvalue = 3, offvalue = 0,bg='#142841',command = test2).place(x=223, y=540)
#d4=Checkbutton(root, text="4",variable=var,onvalue = 4, offvalue = 0,bg='#142841',command = test2).place(x=268, y=540)
#d5=Checkbutton(root, text="5",variable=var,onvalue = 5, offvalue = 0,bg='#142841',command = test2).place(x=313, y=540)
#d6=Checkbutton(root, text="6",variable=var,onvalue = 6, offvalue = 0,bg='#142841',command = test2).place(x=359, y=540)
#d7=Checkbutton(root, text="7",variable=var,onvalue = 7, offvalue = 0,bg='#142841',command = test2).place(x=410, y=540)


#-----------------------------------------------------------------------------------------------------------------------------------------------
#                                                                  BUTTONS
#-----------------------------------------------------------------------------------------------------------------------------------------------

def on_restart():
    root.destroy()
    root2.destroy()
    root4.destroy()
    
btn1 = Button(root, text="Find File", width=44,height=1,bg='#57CAC8', fg='WHITE',font=('ariel 11 bold'), relief=GROOVE, command=model_selection)
btn1.place(x=50, y=80)
btn3 = Button(root, text="Exit", width=8, bg='#142841', fg='white', font=('ariel 11 bold'), relief=GROOVE, command=on_restart)
btn3.place(x=915, y=595)
root4.mainloop()
root2.mainloop()
root.mainloop()




