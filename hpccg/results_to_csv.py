#!/usr/bin/python3
#
# Script to parse the yaml output files of the HPCCG benchmark and
# generate pretty plots.
# by Felix Moessbauer (felix.moessbauer[at]campus.lmu.de)
# 
# usage: python3 plot_results.py
#

import glob
import re
import yaml
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

def import_data (mpi_types):
  results = pd.DataFrame()
  for mpi_type in mpi_types:
    for filename in glob.iglob('results/{}/**/*.yaml'.format(mpi_type), recursive=True):
      print("process", filename)
      parts = re.search('.*\d+/src_(.*)/[^/]+', filename)
      impl = parts.group(1)
      fh = open(filename, 'r')
      data = fh.read()
      ydat = yaml.load(data)
      dffile = pd.io.json.json_normalize(ydat,record_prefix=".")
      dffile['mpi']  = mpi_type
      dffile['impl'] = impl;
      results = pd.concat([results, dffile])
  return results

# Name folders to analyze
mpi = ['cray', 'cray-knl']

data = import_data(mpi);
data.to_csv('data_raw.csv', header=True)
