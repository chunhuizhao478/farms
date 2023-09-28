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

full_path = "./TPV103_Nx/TPV103_Nx720_s2.00_tf0.10_npc5-DataFiles" #100m
# full_path = "./TPV103_Nx/TPV103_Nx1440_s2.00_tf0.10_npc20-DataFiles" #50m
# full_path = "./TPV103_Nx/TPV103_Nx2880_s2.00_tf0.10_npc5-DataFiles" #25m

veldata = read_full_data(full_path+"/top_velo_0.out", 720, 119)

dispdata = read_full_data(full_path+"/top_disp_0.out", 720, 119)

print(np.shape(veldata))

np.savetxt(full_path+'/vel0.txt',veldata)

np.savetxt(full_path+'/disp0.txt',dispdata)