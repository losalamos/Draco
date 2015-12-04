#!/usr/bin/env python

# plot_ipcress_opacity.py
# This program plots opacity from IPCRESS files. It could easily be extended
# to plot other fields from IPCRESS files as well.
# by Alex Long 12/15/2014


# import block
################################################################################
import re

import matplotlib
matplotlib.use('TkAgg')

from numpy import arange, sin, pi, min, max
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

from matplotlib.figure import Figure

import sys 
if sys.version_info[0] < 3:
    from Tkinter import *
else:
    from tkinter import *

from struct import *
import numpy as np
from math import *
################################################################################

# These are the functions that are used to read data from the
# binary IPCRESS file. It also contains a function for interpolating in
# density and temperature. The data locations are specified in 
# cdi_ipcress/doc/IPCRESS_File_Format.pdf

################################################################################
def get_data_for_id(filename, data_start_index, num_entries):
  temp_grid = []
  # "rb" is read binary mode
  with open(filename, "rb") as f:
    f.seek(data_start_index*8)
    for i in range(num_entries):
      word = f.read(8)
      temp_grid.append(unpack('>d', word)[0])
  return temp_grid 
################################################################################


################################################################################
def get_interpolated_data(T_grid, rho_grid, hnu_grid, op_table, target_rho, target_T, print_str):
  n_rho = len(rho_grid)
  n_T = len(T_grid)
  n_hnu = len(hnu_grid)

  # don't allow extrapolation
  if (target_rho  < np.min(rho_grid)):  target_rho = np.min(rho_grid)
  if (target_rho  > np.max(rho_grid)):  target_rho = np.max(rho_grid)
  if (target_T    < np.min(T_grid)):    target_T = np.min(T_grid)
  if (target_T    > np.max(T_grid)):    target_T = np.max(T_grid)
  print("{0} , Target rho: {1} , target T: {2}".format(print_str, target_rho, target_T))

  rho_L = 1000; rho_G =0
  for rho_i, rho in enumerate(rho_grid[:-1]):
    if ( target_rho >= rho and target_rho<=rho_grid[rho_i+1]):
      rho_L = rho_i
      rho_G = rho_i+1
      break

  for T_i, T in enumerate(T_grid[:-1]):
    if ( target_T >= T and target_T<=T_grid[T_i+1]):
      T_L = T_i
      T_G = T_i+1
      break

  #print("Temperature interpolation bounds: {0} {1}".format(T_grid[T_L], T_grid[T_G]))
  #print("Density interpolation bounds: {0} {1}".format(rho_grid[rho_L], rho_grid[rho_G]))

  #get the adjacent rows of the opacity index
  rho_L_T_L = op_table[n_rho*T_L*(n_hnu-1) + rho_L*(n_hnu-1) : n_rho*T_L*(n_hnu-1) + rho_L*(n_hnu-1) + (n_hnu-1) ]
  rho_L_T_G = op_table[n_rho*T_G*(n_hnu-1) + rho_L*(n_hnu-1) : n_rho*T_G*(n_hnu-1) + rho_L*(n_hnu-1) + (n_hnu-1) ]
  rho_G_T_L = op_table[n_rho*T_L*(n_hnu-1) + rho_G*(n_hnu-1) : n_rho*T_L*(n_hnu-1) + rho_G*(n_hnu-1) + (n_hnu-1) ]
  rho_G_T_G = op_table[n_rho*T_G*(n_hnu-1) + rho_G*(n_hnu-1) : n_rho*T_G*(n_hnu-1) + rho_G*(n_hnu-1) + (n_hnu-1) ]

  interp_op = []
  #interpolate for each frequency point
  for i in range(n_hnu-1):
    #if (rho_L_T_L[i] < 1.0e-10) or (rho_L_T_G[i] < 1.0e-10) or (rho_G_T_L[i] < 1.0e-10) or (rho_G_T_G[i] < 1.0e10):
    #  interp_op.append(1.0e-10)
    #print("{0} {1} {2} {3}" .format(rho_L_T_L[i], rho_L_T_G[i], rho_G_T_L[i], rho_G_T_G[i]))
    log_op_T_L  = log(rho_L_T_L[i]) + log(target_rho/rho_grid[rho_L]) / log(rho_grid[rho_G]/rho_grid[rho_L]) * log(rho_G_T_L[i]/rho_L_T_L[i])
    log_op_T_G  = log(rho_L_T_G[i]) + log(target_rho/rho_grid[rho_L]) / log(rho_grid[rho_G]/rho_grid[rho_L]) * log(rho_G_T_G[i]/rho_L_T_G[i])
    log_op =  log_op_T_L + log(target_T/T_grid[T_L]) / log(T_grid[T_G]/T_grid[T_L]) * (log_op_T_G - log_op_T_L)
    interp_op.append(exp(log_op))

  print("hnu(keV)      opacity(sq_cm/g)     opacity(1/cm)")
  for i, hnu in enumerate(hnu_grid[:-1]):
    print("{0}   {1}   {2}".format( 0.5*(hnu + hnu_grid[i+1]), interp_op[i], interp_op[i]*target_rho))
  return interp_op 
###############################################################################

###############################################################################
def read_information_from_file(ipcress_file):

  word_array = []
  with open(ipcress_file, "rb") as f:
    for i in range(26):
      word = f.read(8)
      if not word:
          break
      word_array.append(word)
      #print(int(unpack('>d', word)[0]))

  title = word_array[0]
  toc_int= []
  offset = 2
  for i in range(offset,offset+24):
    toc_int.append( int(unpack('>d', word_array[i])[0]))

  n_data_records = toc_int[14]
  mxrec = toc_int[1] - toc_int[0]
  mxkey = toc_int[16]
  #print("Number of data records: {0}".format(n_data_records))
  #print("Beginnging of data: {0}".format(toc_int[0]))
  #print("Max records: {0} , max search keys: {1}".format(mxrec, mxkey))

  property = []

  ds = []
  dfo = []
  tdf = []
  num_mats = 0
  mat_ids= []
  with open(ipcress_file, "rb") as f:
    # Read in array that lists the data sizes in this file
    f.seek(toc_int[0]*8)
    #print("Table of data sizes")
    for i in range(n_data_records):
      word = f.read(8)
      ds.append(int(unpack('>d', word)[0]))

    # Read in array gives the offsets between data 
    f.seek(toc_int[1]*8)
    #print("Table of data file offesets")
    for i in range(n_data_records):
      word = f.read(8)
      dfo.append(int(unpack('>d', word)[0]))

    # Read in material IDs present in this file
    f.seek(dfo[0]*8)
    #print("Table of  material identifiers")
    word = f.read(8)
    num_mats = int(unpack('>d', word)[0])
    for i in range(num_mats):
      word = f.read(8)
      mat_ids.append( int(unpack('>d', word)[0]))

    # Read in list of properties in this file available for each material
    # entries in this table are 24 bytes each
    f.seek(toc_int[10]*8)
    #print("Table of data fields for each  material")
    word = f.read(72) #ignore the first 72 bytes, they don't contain useful information
    for i in range(1,toc_int[14]):
      #oredering is "matID" "data type" "fill"
      temp_property = []
      for j in range(mxkey):
        three_string = []
        three_string.append( f.read(8))
        three_string.append( f.read(8))
        three_string.append( f.read(8))
        if (j==0): temp_property.append(three_string[2].strip() )
        elif (j==1): temp_property.append(three_string[0].strip())
        else: temp_property.append(i) #index of data table containing values
      property.append(temp_property)

  materials = []
  for m in range(num_mats):
    materials.append([ m, mat_ids[m]])

  print("{0} materials in file".format(num_mats))
  for i in range(num_mats):
    print("  Matieral ID: {0}".format(mat_ids[i]))

  print("List of available properties")
  for i in property:
    print(i)
  
  #return the list of available properties, data file offsets and data sizes
  return materials, property, dfo, ds
################################################################################

################################################################################
# Checks to see if there are any zeros in the opcaity data--zero data is 
# difficult to handle and for now we are going to ignore data sets that contain
# zeros and print an error message
def check_valid_data(opacity_grid):
  for item in opacity_grid:
    if (item != 0.0):
      return True
  return False
################################################################################



###############################################################################
#  MAIN PROGRAM
###############################################################################


# make sure an IPCRESS file is specified
if (len(sys.argv) != 2):
  print("Useage: {0} <path to ipcress file>".format(sys.argv[0]))
  sys.exit()

ipcress_file = sys.argv[1]

# regular expressions for properties
re_ramg = re.compile("ramg")
re_pmg = re.compile("pmg")
re_rsmg = re.compile("rsmg")
re_rtmg = re.compile("rtmg")
re_tgrid = re.compile("tgrid")
re_hnu = re.compile("hnugrid")
re_rho = re.compile("rgrid")


# set up TK plot
root = Tk()
root.wm_title("Data from IPCRESS File {0}".format(ipcress_file))

#load data from IPCRESS file
# dfo is the array of data file offsets, ds is the array of data sizes
materials, property_list, dfo, ds = read_information_from_file(ipcress_file)

#build dictionary of data, keys are "property_matID"
table_key_dict = {}
for mat in materials:
  mat_ID = mat[1]
  for prop_i, prop in enumerate(property_list):
    if re_ramg.search(prop[1])  and int(prop[0]) == mat_ID:
      table_key_dict["{0}_{1}".format("ramg", mat_ID)] = get_data_for_id( ipcress_file, dfo[prop_i+1], ds[prop_i+1])
    if re_rsmg.search(prop[1])  and int(prop[0]) == mat_ID:
      table_key_dict["{0}_{1}".format("rsmg", mat_ID)] = get_data_for_id( ipcress_file, dfo[prop_i+1], ds[prop_i+1])
    if re_pmg.search(prop[1])  and int(prop[0]) == mat_ID:
      table_key_dict["{0}_{1}".format("pmg", mat_ID)] = get_data_for_id( ipcress_file, dfo[prop_i+1], ds[prop_i+1])
    if re_tgrid.search(prop[1])  and int(prop[0]) == mat_ID:
      table_key_dict["{0}_{1}".format("tgrid", mat_ID)] = get_data_for_id( ipcress_file, dfo[prop_i+1], ds[prop_i+1])
    if re_hnu.search(prop[1])  and int(prop[0]) == mat_ID:
      table_key_dict["{0}_{1}".format("hnugrid", mat_ID)] = get_data_for_id( ipcress_file, dfo[prop_i+1], ds[prop_i+1])
    if re_rho.search(prop[1])  and int(prop[0]) == mat_ID:
      table_key_dict["{0}_{1}".format("rgrid", mat_ID)] = get_data_for_id( ipcress_file, dfo[prop_i+1], ds[prop_i+1])


# select the first material in the file as default and set as a TK variable
selected_ID = IntVar()  # make this a TK variable
selected_ID.set(materials[0][1]) # set TK variable

# add figuree
f = Figure(figsize=(5,4), dpi=100)
a = f.add_subplot(111)

# set values for the first material (selected)
T_grid = table_key_dict["{0}_{1}".format("tgrid", selected_ID.get())]
rho_grid = table_key_dict["{0}_{1}".format("rgrid", selected_ID.get())]
hnu_grid = table_key_dict["{0}_{1}".format("hnugrid", selected_ID.get())]
mgr_grid = table_key_dict["{0}_{1}".format("ramg", selected_ID.get())]
mgp_grid = table_key_dict["{0}_{1}".format("pmg", selected_ID.get())]
mgs_grid = table_key_dict["{0}_{1}".format("rsmg", selected_ID.get())]

#get log average rho and T for initial plot
avg_rho = exp(0.5*(log(max(rho_grid)) + log(min(rho_grid))))
avg_rho = round(avg_rho, -int(floor(log10(avg_rho))-1))
avg_T = exp(0.5*(log(max(T_grid)) + log(min(T_grid))))
avg_T = round(avg_T, -int(floor(log10(avg_T))-1)) 

# setup and initialize Tk values that can be tied to entry fields
# set temperature and density to log average value
target_rho = DoubleVar()
target_rho.set(avg_rho)
target_T = DoubleVar()
target_T.set(avg_T)
op_type = IntVar()
op_type.set(1)

# get interpolated opacity data for this temperature and density
# first check to see if the data is valid (non-zero)
mgr_valid = check_valid_data(mgr_grid)
mgp_valid = check_valid_data(mgp_grid)
mgs_valid = check_valid_data(mgs_grid)

# if valid, interpolate data at target rho and target T
if (mgr_valid):
  mgr_interp = get_interpolated_data(T_grid, rho_grid, hnu_grid, mgr_grid, \
                target_rho.get(), target_T.get(), "multigroup absorption Rosseland")
if (mgp_valid):
  mgp_interp = get_interpolated_data(T_grid, rho_grid, hnu_grid, mgp_grid, \
                target_rho.get(), target_T.get(), "multigroup absorption Planckian")
if (mgs_valid):
  mgs_interp = get_interpolated_data(T_grid, rho_grid, hnu_grid, mgs_grid, \
                target_rho.get(), target_T.get(), "multigroup scattering")

# plotting data arrays
opr_data = []  # Rosseland absorption
opp_data = []  # Planck absorption
ops_data = []  # Rosseland scattering
hnu_data = []  # hnu plot data

# add points to plotting data twice to get the bar structure and find max and min opacity
min_opacity = 1.0e10
max_opacity = 0.0
for hnu_i, hnu in enumerate(hnu_grid[:-1]):
  hnu_data.append(hnu)
  hnu_data.append(hnu_grid[hnu_i+1])
  if (mgr_valid):
    opr_data.append(mgr_interp[hnu_i])
    opr_data.append(mgr_interp[hnu_i])
    min_opacity = min([min(opr_data), min_opacity])
    max_opacity = max([max(opr_data), max_opacity])

  if (mgp_valid):
    opp_data.append(mgp_interp[hnu_i])
    opp_data.append(mgp_interp[hnu_i])
    min_opacity = min([min(opp_data), min_opacity])
    max_opacity = max([max(opp_data), max_opacity])

  if (mgs_valid):
    ops_data.append(mgs_interp[hnu_i])
    ops_data.append(mgs_interp[hnu_i])
    min_opacity = min([min(ops_data), min_opacity])
    max_opacity = max([max(ops_data), max_opacity])

# Plot data with the initial settings
a.set_xscale('log') 
a.set_yscale('log') 
a.set_xlabel('hnu (keV)')
a.set_ylabel('opacity (sq. cm/g)')
if (mgr_valid):
  a.plot(hnu_data, opr_data, 'b-', label = "{0} Rosseland Absorption".format(selected_ID.get()))
else:
  print("ERROR: Invalid multigroup Rosseland absorption data (all zeros)") 

if (mgs_valid):
  a.plot(hnu_data, ops_data, 'r-', label = "{0} Rosseland Scattering".format(selected_ID.get()))
else:
  print("ERROR: Invalid multigroup scattering data (all zeros)")

a.set_xlim([0.9*min(hnu_grid), 1.1*max(hnu_grid)]) 
a.set_ylim([ 0.7*min_opacity, 1.3*max_opacity]) 

a.legend(loc='best')

canvas = FigureCanvasTkAgg(f, master=root)
canvas.show()
canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

toolbar = NavigationToolbar2TkAgg( canvas, root )
toolbar.update()
canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)


###############################################################################
# Plot loop executed when the "Plot" button is pressed
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
  mgr_grid = table_key_dict["{0}_{1}".format("ramg", selected_ID.get())]
  mgp_grid = table_key_dict["{0}_{1}".format("pmg", selected_ID.get())]
  # scattering is only availabe with Rosseland weighting
  mgs_grid = table_key_dict["{0}_{1}".format("rsmg", selected_ID.get())]
  
  # get interpolated opacity data for this temperature and density
  # first check to see if the data is valid (non-zero)
  mgr_valid = check_valid_data(mgr_grid)
  mgp_valid = check_valid_data(mgp_grid)
  mgs_valid = check_valid_data(mgs_grid)

  # if valid, interpolate data at target rho and target T
  if (mgr_valid):
    mgr_interp = get_interpolated_data(T_grid, rho_grid, hnu_grid, mgr_grid, \
                  target_rho.get(), target_T.get(), "multigroup absorption Rosseland")
  if (mgp_valid):
    mgp_interp = get_interpolated_data(T_grid, rho_grid, hnu_grid, mgp_grid, \
                  target_rho.get(), target_T.get(), "multigroup absorption Planckian")
  if (mgs_valid):
    mgs_interp = get_interpolated_data(T_grid, rho_grid, hnu_grid, mgs_grid, \
                  target_rho.get(), target_T.get(), "multigroup scattering")

  # plotting data arrays
  opr_data = []  # Rosseland absorption
  opp_data = []  # Planck absorption
  ops_data = []  # Rosseland scattering
  hnu_data = []  # hnu plot data

  # add points to plotting data twice to get the bar structure and find max and min opacity
  min_opacity = 1.0e10
  max_opacity = 0.0
  for hnu_i, hnu in enumerate(hnu_grid[:-1]):
    hnu_data.append(hnu)
    hnu_data.append(hnu_grid[hnu_i+1])

    if (mgr_valid):
      opr_data.append(mgr_interp[hnu_i])
      opr_data.append(mgr_interp[hnu_i])
      min_opacity = min([min(opr_data), min_opacity])
      max_opacity = max([max(opr_data), max_opacity])

    if (mgp_valid):
      opp_data.append(mgp_interp[hnu_i])
      opp_data.append(mgp_interp[hnu_i])
      min_opacity = min([min(opp_data), min_opacity])
      max_opacity = max([max(opp_data), max_opacity])

    if (mgs_valid):
      ops_data.append(mgs_interp[hnu_i])
      ops_data.append(mgs_interp[hnu_i])
      min_opacity = min([min(ops_data), min_opacity])
      max_opacity = max([max(ops_data), max_opacity])
  
  a.set_xscale('log') 
  a.set_yscale('log') 
  a.set_xlabel('hnu (keV)')
  a.set_ylabel('opacity (sq. cm/g)')
  # use label for Planck or Rosseland
  if (op_type.get() == 1 and mgr_valid ):
    a.plot(hnu_data, opr_data, 'b-',  label = "{0} Rosseland Absorption".format(selected_ID.get()))
  elif(op_type.get() == 2 and mgp_valid):
    a.plot(hnu_data, opp_data, 'g-', label = "{0} Planckian Absorption".format(selected_ID.get()))
  elif(op_type.get() == 3 and mgr_valid and mgp_valid):
    a.plot(hnu_data, opr_data, 'b-',  label = "{0} Rosseland Absorption".format(selected_ID.get()))
    a.plot(hnu_data, opp_data, 'g-', label = "{0} Planckian Absorption".format(selected_ID.get()))
  else:
    if (not mgr_valid):
      print("ERROR: Invalid multigroup Rosseland absorption data (all zeros)")
    if (not mgp_valid):
      print("ERROR: Invalid multigroup Planckian absorption data (all zeros)")
  if (mgs_valid):
    a.plot(hnu_data, ops_data, 'r-', label = "{0} Rosseland Scattering".format(selected_ID.get()))
  else:
    print("ERROR: Invalid multigroup scattering data (all zeros)")

  a.legend(loc='best')
  a.set_xlim([0.9*min(hnu_grid), 1.1*max(hnu_grid)]) 
  a.set_ylim([ 0.7*min_opacity, 1.3*max_opacity]) 
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


L_opacity = Label(root, text = "Opacity Weighting Method")
L_opacity.pack(anchor = W)
Radiobutton(root, text="Rosseland", padx = 20, variable=op_type, value=1).pack(anchor=W)
Radiobutton(root, text="Planckian",  padx = 20, variable=op_type, value=2).pack(anchor=W)
Radiobutton(root, text="Show Both",  padx = 20, variable=op_type, value=3).pack(anchor=W)


button = Button(master=root, text='Quit', command=_quit)
button.pack(side=BOTTOM)


button = Button(root, text='Plot', command=plot_op)
button.pack( side=BOTTOM ) 

mainloop()
