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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

def import_data (mpi_type):
  results = defaultdict(dict)
  for filename in glob.iglob('results/{}/**/*.yaml'.format(mpi_type), recursive=True):
    parts = re.search('.*\d+/src_(.*)/[^/]+', filename)
    impl = parts.group(1)
    fh = open(filename, 'r')
    data = fh.read()
    ydat = yaml.load(data)
    ranks = ydat['Parallelism']['Number of MPI ranks']
    if not impl in results:
      results[impl] = defaultdict(list)
    results[impl][ranks].append(ydat['MFLOPS Summary']['Total'])
  return results

def flatten_results (results):
  flat_results = defaultdict(dict) 
  for item in results.items(): # mpi versions
    mpi, data = item
    for impl, res in data.items(): # impls
      for node in res.items():     # number of nodes
        units, mflops = node
        mean = np.mean(mflops)
        if not impl in flat_results[mpi]:
          flat_results[mpi][impl] = defaultdict(dict)
        flat_results[mpi][impl][units] = mean
  return(flat_results)

# Perform analysis
mpi_impls = ['intel', 'ibm']
results = {}

for mpi_impl in mpi_impls:
  results[mpi_impl] = import_data(mpi_impl)

flat_results = flatten_results(results)

# melt data for plotting
df = pd.Panel(flat_results).to_frame().reset_index()

# Plotting results
for mpi_impl in mpi_impls:
  for (name, group) in df.groupby("minor"):
    plt.scatter(x=group['major'], y=group[mpi_impl], label=name)
  
  #plt.yscale("log")
  plt.xscale("log")
  ticks = [1,2,4,8,16,28] 
  plt.xticks(ticks, ticks)
  plt.xlabel("units")
  plt.ylabel("MFLOPS")
  plt.title("HPCCG Performance using {} MPI".format(mpi_impl))
  plt.legend()
  plt.grid();
  plt.savefig("results/hpccg_performance_{}_mpi.svg".format(mpi_impl), format="svg")
#  plt.show()

