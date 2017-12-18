#!/usr/bin/env python

import sys
import csv
import plotly
import plotly.plotly
import plotly.figure_factory
import numpy
import pandas

y = []
gr = []
with open(sys.argv[1] if len(sys.argv) > 1 else 'trace.csv', 'rb') as f:
    for row in csv.reader(f, delimiter=';'):
        if row[1] not in ('smooth_inner', 'smooth_res_inner'):
            continue
        y  += [float(row[5]) / float(row[3]) / float(1024) / float(1024) / float(1024)]
        gr += [row[1]]

frames = pandas.DataFrame(dict(Flops = y, Phase = gr))
fig = plotly.figure_factory.create_violin(frames, data_header='Flops', group_header='Phase', height=786, width=1024)
plotly.offline.plot(fig, filename='Flops.html', auto_open=False)
