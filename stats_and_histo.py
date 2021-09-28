#!/usr/bin/env python3
import random as rnd
import numpy as np
import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt
import graph_tool as gt
from graph_tool.all import *
import pandas as pd
import time
import scipy.special
import scipy.stats


#The script contains a function that saves arrays as .npz files

def out_to_tex(path, array, unshuffled, name):
    np.savez('{0}/{1}.npz'.format(path,name), array=array, unshuffled=unshuffled)
    

