#!/usr/bin/env python

# plot_ipcress_opacity.py
# This program plots opacity from IPCRESS files. It could easily be extended
# to plot other fields from IPCRESS files as well.
# NOTE: Currently, this will not run on CCS machines because they do not have
# backend_tkagg installed 

# by Alex Long 12/15/2014

import re

import matplotlib
matplotlib.use('TkAgg')

from numpy import arange, sin, pi, min, max
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

from matplotlib.figure import Figure
import ipcress_reader_functions as op_read

import sys 
if sys.version_info[0] < 3:
    from Tkinter import *
else:
    from tkinter import *


# make sure an IPCRESS file is specified
if (len(sys.argv) != 2):
  print("Useage: {0} <path to ipcress file>".format(sys.argv[0]))
  sys.exit()

ipcress_file = sys.argv[1]

# regular expressions for properties
re_ramg = re.compile("ramg")
re_rsmg = re.compile("rsmg")
re_rtmg = re.compile("rtmg")
re_tgrid = re.compile("tgrid")
re_hnu = re.compile("hnugrid")
re_rho = re.compile("rgrid")


root = Tk()
root.wm_title("Data from IPCRESS File {0}".format(ipcress_file))

#load data from IPCRESS file
# dfo is the array of data file offsets, ds is the array of data sizes
materials, property_list, dfo, ds = op_read.read_information_from_file(ipcress_file)

#build dictionary of data
table_key_dict = {}
for mat in materials:
  mat_ID = mat[1]
  for prop_i, prop in enumerate(property_list):
    if re_ramg.search(prop[1])  and int(prop[0]) == mat_ID:
      table_key_dict["{0}_{1}".format("ramg", mat_ID)] = op_read.get_data_for_id( ipcress_file, dfo[prop_i+1], ds[prop_i+1])
    if re_rsmg.search(prop[1])  and int(prop[0]) == mat_ID:
      table_key_dict["{0}_{1}".format("rsmg", mat_ID)] = op_read.get_data_for_id( ipcress_file, dfo[prop_i+1], ds[prop_i+1])
    if re_tgrid.search(prop[1])  and int(prop[0]) == mat_ID:
      table_key_dict["{0}_{1}".format("tgrid", mat_ID)] = op_read.get_data_for_id( ipcress_file, dfo[prop_i+1], ds[prop_i+1])
    if re_hnu.search(prop[1])  and int(prop[0]) == mat_ID:
      table_key_dict["{0}_{1}".format("hnugrid", mat_ID)] = op_read.get_data_for_id( ipcress_file, dfo[prop_i+1], ds[prop_i+1])
    if re_rho.search(prop[1])  and int(prop[0]) == mat_ID:
      table_key_dict["{0}_{1}".format("rgrid", mat_ID)] = op_read.get_data_for_id( ipcress_file, dfo[prop_i+1], ds[prop_i+1])
     
# set Tk values that can be tied to entry fields
target_rho = DoubleVar()
target_rho.set(1.0)
target_T = DoubleVar()
target_T.set(1.0)
# select the first material in the file as default
selected_ID = IntVar()
selected_ID.set(materials[0][1])


f = Figure(figsize=(5,4), dpi=100)
a = f.add_subplot(111)

# set values for the first material (selected)
T_grid = table_key_dict["{0}_{1}".format("tgrid", selected_ID.get())] 
rho_grid = table_key_dict["{0}_{1}".format("rgrid", selected_ID.get())] 
hnu_grid = table_key_dict["{0}_{1}".format("hnugrid", selected_ID.get())] 
mg_grid = table_key_dict["{0}_{1}".format("ramg", selected_ID.get())]
mgs_grid = table_key_dict["{0}_{1}".format("rsmg", selected_ID.get())]

# get interpolated opacity data for this temperature and density
mg_interp = op_read.get_interpolated_data(T_grid, rho_grid, hnu_grid, mg_grid, target_rho.get(), target_T.get())
mgs_interp = op_read.get_interpolated_data(T_grid, rho_grid, hnu_grid, mgs_grid, target_rho.get(), target_T.get())
op_data = []
ops_data = []
hnu_data = []
for hnu_i, hnu in enumerate(hnu_grid[:-1]):
  hnu_data.append(hnu)
  hnu_data.append(hnu_grid[hnu_i+1])
  op_data.append(mg_interp[hnu_i])
  op_data.append(mg_interp[hnu_i])
  ops_data.append(mgs_interp[hnu_i])
  ops_data.append(mgs_interp[hnu_i])

a.set_xscale('log') 
a.set_yscale('log') 
a.set_xlabel('hnu (keV)')
a.set_ylabel('opacity (sq. cm/g)')
a.plot(hnu_data, op_data, label = "{0} Rosseland Absorption".format(selected_ID.get()))
a.plot(hnu_data, ops_data, label = "{0} Rosseland Scattering".format(selected_ID.get()))
a.set_xlim([0.9*min(hnu_grid), 1.1*max(hnu_grid)]) 
a.set_ylim([ 0.7*min( [min(op_data), min(ops_data)]), 1.3*max( [max(op_data), max(ops_data)])]) 
a.legend(loc='best')

canvas = FigureCanvasTkAgg(f, master=root)
canvas.show()
canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

toolbar = NavigationToolbar2TkAgg( canvas, root )
toolbar.update()
canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
 

###############################################################################
def plot_op():
  f.clf()
  
  a = f.add_subplot(111)

  #make sure selected material is available, if not found set to first material
  mat_found = False
  for mat in materials:
    if selected_ID.get() == mat[1]:
      mat_found=True
  if mat_found == False: selected_ID.set(materials[0][1])

  # get values associated with this ID from dictionary
  T_grid = table_key_dict["{0}_{1}".format("tgrid", selected_ID.get())] 
  rho_grid = table_key_dict["{0}_{1}".format("rgrid", selected_ID.get())] 
  hnu_grid = table_key_dict["{0}_{1}".format("hnugrid", selected_ID.get())] 
  mg_grid = table_key_dict["{0}_{1}".format("ramg", selected_ID.get())]
  mgs_grid = table_key_dict["{0}_{1}".format("rsmg", selected_ID.get())]

  
  # get interpolated opacity data for this temperature and density
  mg_interp = op_read.get_interpolated_data(T_grid, rho_grid, hnu_grid, mg_grid, target_rho.get(), target_T.get())
  mgs_interp = op_read.get_interpolated_data(T_grid, rho_grid, hnu_grid, mgs_grid, target_rho.get(), target_T.get())
  # make the piecewise constant data to be plotted
  op_data = []
  ops_data = []
  hnu_data = []
  for hnu_i, hnu in enumerate(hnu_grid[:-1]):
    hnu_data.append(hnu)
    hnu_data.append(hnu_grid[hnu_i+1])
    op_data.append(mg_interp[hnu_i])
    op_data.append(mg_interp[hnu_i])
    ops_data.append(mgs_interp[hnu_i])
    ops_data.append(mgs_interp[hnu_i])
  
  a.set_xscale('log') 
  a.set_yscale('log') 
  a.set_xlabel('hnu (keV)')
  a.set_ylabel('opacity (sq. cm/g)')
  a.plot(hnu_data, op_data, label = "{0} Rosseland Absorption".format(selected_ID.get()))
  a.plot(hnu_data, ops_data, label = "{0} Rosseland Scattering".format(selected_ID.get()))
  a.legend(loc='best')
  a.set_xlim([0.9*min(hnu_grid), 1.1*max(hnu_grid)]) 
  a.set_ylim([ 0.7*min( [min(op_data), min(ops_data)]), 1.3*max( [max(op_data), max(ops_data)])]) 
  canvas.show()
  canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
###############################################################################


# set up entry fields and their associated labels

L_rho = Label(root, text = "Density ({0} g/cc - {1} g/cc)".format(min(rho_grid), max(rho_grid)))
L_rho.pack(anchor = W)
target_rho_field = Entry(root, textvariable=target_rho)
target_rho_field.pack( anchor = W ) 

L_T = Label(root, text = "Temperature ({0} keV - {1} keV)".format(min(T_grid), max(T_grid)))
L_T.pack(anchor = W)
target_T_field = Entry(root, textvariable=target_T)
target_T_field.pack( anchor = W ) 

material_str = ""
for mat in materials: material_str = "{0} {1}".format(material_str, mat[1])
L_material =Label(root, text = "Materials in file: {0} ".format(material_str))
L_material.pack(anchor = W)
material_field = Entry(root, textvariable=selected_ID)
material_field.pack( anchor = W ) 


###############################################################################
def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate
###############################################################################

button = Button(master=root, text='Quit', command=_quit)
button.pack(side=BOTTOM)


button = Button(root, text='Plot', command=plot_op)
button.pack( side=BOTTOM ) 

mainloop()
