# ipcress_reader_functions.py
# This program contains the functions that are used to read data from the
# binary IPCRESS file. It also contains a function for interpolating in
# density and temperature. The data locations are specified in 
# cdi_ipcress/doc/IPCRESS_File_Format.pdf
# by Alex Long 12/15/2014

from struct import *
import numpy as np
from math import *


###############################################################################
def get_data_for_id(filename, data_start_index, num_entries):
  temp_grid = []
  with open(filename, "rb") as f:
    f.seek(data_start_index*8)
    for i in range(num_entries):
      word = f.read(8)
      temp_grid.append(unpack('>d', word)[0])
  return temp_grid 
###############################################################################


###############################################################################
def get_interpolated_data(T_grid, rho_grid, hnu_grid, op_table, target_rho, target_T):
  n_rho = len(rho_grid)
  n_T = len(T_grid)
  n_hnu = len(hnu_grid)

  #don't allow extrapolation
  if (target_rho  < np.min(rho_grid)):  target_rho = np.min(rho_grid)
  if (target_rho  > np.max(rho_grid)):  target_rho = np.max(rho_grid)
  if (target_T    < np.min(T_grid)):    target_T = np.min(T_grid)
  if (target_T    > np.max(T_grid)):    target_T = np.max(T_grid)
  print("Target rho: {0} , target T: {1}".format(target_rho, target_T))

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
    log_op_T_L  = log(rho_L_T_L[i]) + log(target_rho/rho_grid[rho_L]) / log(rho_grid[rho_G]/rho_grid[rho_L]) * log(rho_G_T_L[i]/rho_L_T_L[i])
    log_op_T_G  = log(rho_L_T_G[i]) + log(target_rho/rho_grid[rho_L]) / log(rho_grid[rho_G]/rho_grid[rho_L]) * log(rho_G_T_G[i]/rho_L_T_G[i])
    log_op =  log_op_T_L + log(target_T/T_grid[T_L]) / log(T_grid[T_G]/T_grid[T_L]) * (log_op_T_G - log_op_T_L)
    interp_op.append(exp(log_op))
  print(interp_op)
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
###############################################################################
