#!/usr/bin/env python
#-----------------------------*-python-*----------------------------------------#
# file   src/cdi_ipcress/python/demo_read_data.py
# author Alex Long <along@lanl.gov>
# date   Friday, October 28, 2016, 5:44 pm
# brief  A script that uses the ipcress_reader.py functions to get gray opacity
#        data for a range of temperatures.
# note   Copyright (C) 2016, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
import sys
import os
import ipcress_reader as ip_reader

# make sure an IPCRESS file is specified
if (len(sys.argv) != 2):
  print("Usage: {0} <path to ipcress file>".format(sys.argv[0]))
  sys.exit()

# The property names in the output map as:
# rgray = Gray Rosseland total
# ragray = Gray Rosseland absorption
# ramg = MG Rosseland absorption
# p = Plank total? (not used in draco)
# pgray = Gray Planck absorption

ipcress_file = sys.argv[1]

# make sure ipcress file exists
if (not os.path.exists(ipcress_file)):
  print("File at {0} does not exist... exiting".format(ipcress_file))
  sys.exit()

# get data dictionary and material IDs from file
ipcress_data, materials = \
  ip_reader.get_property_map_from_ipcress_file(ipcress_file)

# print the available data in this file
print("Available data sets in {0}".format(ipcress_file))
print(ipcress_data.keys())
print("\n")

# check for data set I want to use
mat_ID = 10001
if "{0}_{1}".format("rgray",  mat_ID) in ipcress_data:
  print("File has Rosseland gray data for {0}".format(mat_ID))
else:
  print("File does not have Rosseland gray data for {0}".format(mat_ID))
  sys.exit(0)
print("\n")

# get Rosseland gray data table
rgray_data = ipcress_data["{0}_{1}".format("rgray", mat_ID)]

# get data need for interpolating gray Rosseland data
rho_grid = ipcress_data["{0}_{1}".format("rgrid",  mat_ID)]
T_grid = ipcress_data["{0}_{1}".format("tgrid",  mat_ID)]

# interpolate for a list of material temperatures at a given density
T_list = [0.05, 0.1, 0.5, 0.6, 1.25, 1.461, 4.0] #keV
target_rho = 1.313 # g/cc
op_interp = []
for target_T in T_list:
  op_interp.append(ip_reader.interpolate_gray_opacity_data(T_grid, rho_grid, \
    rgray_data, target_rho, target_T))
print("\n")

print("Interpolated gray Rosseland opacities at rho={0} g/cc".format(target_rho))
print("T (kev)    sigma (cm^2/g)   sigma (1/cm)")
for i_T, T in enumerate(T_list):
  print("{0:0.2e}   {1:0.6e}     {2:0.6e}".format(T, op_interp[i_T], \
    op_interp[i_T]*target_rho))

'''
# Multigroup example

hnu_grid = ipcress_data["{0}_{1}".format("hnugrid",  mat_ID)]
pmg_data =  ipcress_data["{0}_{1}".format("pmg",  mat_ID)]
interp_mgp = ip_reader.interpolate_mg_opacity_data(T_grid, rho_grid, hnu_grid, pmg_data, \
target_rho, 1.25)
print(interp_mgp)
'''
