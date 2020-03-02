#!/usr/bin/env python
#-----------------------------*-python-*----------------------------------------#
# file   src/cdi_ipcress/python/demo_modify_ipcress_file.py
# author Alex Long <along@lanl.gov>
# date   Monday, March 18, 2019, 5:44 pm
# brief  A script that uses the ipcress_reader.py functions to write new opcaity
#        values to an IPCRESS file
#        data for a range of temperatures.
# note   Copyright (C) 2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
import sys
import os
import numpy as np
import ipcress_reader as ip_reader
from  shutil import copyfile

# make sure an IPCRESS file is specified
if (len(sys.argv) != 2):
  print("Usage: {0} <path to ipcress file>".format(sys.argv[0]))
  sys.exit()

# The property names in the output map as:
# rgray = Gray Rosseland total
# ragray = Gray Rosseland absorption
# ramg = MG Rosseland absorption
# pgray = Gray Planck absorption

ipcress_file = sys.argv[1]

# make sure ipcress file exists
if (not os.path.exists(ipcress_file)):
  print("File at {0} does not exist... exiting".format(ipcress_file))
  sys.exit()

# get data dictionary and material IDs from file
ipcress_data, materials = \
  ip_reader.get_property_map_from_ipcress_file(ipcress_file)

material_ID = "10001" # must be a string
scat_mat_property = "rsmg" # rosseland scattering
abs_mat_property = "ramg" # absorption scattering
mgsr_grid = ipcress_data["{0}_{1}".format(scat_mat_property, material_ID)]
mgar_grid = ipcress_data["{0}_{1}".format(abs_mat_property, material_ID)]

# get the number of grid points in the opacity data (your new data need to be
# the same size)
n_grid_points = len(mgsr_grid)

# define a new scattering cross section
new_scattering_value = 3.0 # cm**2 / g
new_rsmg = np.full(n_grid_points, new_scattering_value)

# modify the absorption data in place with a limiter
opacity_limit = 10000.0
new_mgar_grid = np.zeros(n_grid_points)
for i, opac in enumerate(mgar_grid):
  if opac > opacity_limit:
    opac = opacity_limit
  new_mgar_grid[i] = opac

# copy the ipcress file so you can modify it
new_ipcress_file = "{0}_new".format(ipcress_file)
copyfile(ipcress_file, new_ipcress_file)

# write the new scattering data to the file
ip_reader.write_information_to_file(new_ipcress_file, \
  material_ID, scat_mat_property, new_rsmg)

# write the limited absorption data to the file
ip_reader.write_information_to_file(new_ipcress_file, \
  material_ID, abs_mat_property, new_mgar_grid)

# get data dictionary and material IDs from the new file
new_ipcress_data, new_materials = \
  ip_reader.get_property_map_from_ipcress_file(new_ipcress_file)

# print the new scattering grid
print("New scattering data")
print(new_ipcress_data["{0}_{1}".format(scat_mat_property, material_ID)])
