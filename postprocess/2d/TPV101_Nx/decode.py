from __future__ import print_function
from __future__ import division
from glob import glob
import sys
import struct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from plot_utils import *

def read_full_data(path, nx, nt):
  # Only reads binary file (float32)
  file = open(path, 'r')
  data=np.fromfile(file,dtype=np.float32)
  print(np.shape(data))
  data=data[:nt*nx]
  data=data.reshape((nt,nx))
  file.close()
  return data

# full_path = "./TPV101_Nx360_s2.00_tf0.35_npc1-DataFiles" #200m
# full_path = "./TPV101_Nx720_s2.00_tf0.35_npc1-DataFiles" #100m
full_path = "./TPV101_Nx1440_s2.00_tf0.35_npc1-DataFiles" #50m

veldata = read_full_data(full_path+"/top_velo_1.out", 1440, 149)

print(np.shape(veldata))

np.savetxt(full_path+'/vel1.txt',veldata)