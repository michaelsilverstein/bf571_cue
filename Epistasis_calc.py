#!/usr/bin/env python3

import sys
import pickle
import numpy as np
import cobra
from cobra.flux_analysis import flux_variability_analysis
import cue_tools as ct

spp = sys.argv[1]
url = 'http://bigg.ucsd.edu/static/models/%s.mat' %spp

model = ct.loadModelURL(url, 'mat')
epi,epi_matrix = ct.computeEpistasis(model, export_matrix = True)

with open('%s_distribution.pickle' %spp,'wb') as f:
    pickle.dump(epi, f)
    
with open('%s_matrix.pickle' %spp,'wb') as f:
    pickle.dump(epi_matrix, f)

