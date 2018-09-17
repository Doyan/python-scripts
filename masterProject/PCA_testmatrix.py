# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 09:42:15 2018

Script for doinng a PCA analysis on the testmatrix we ran during the summer. 

@author: Gabriel
"""

import numpy as np

import pandas as pd

import matplotlib.pyplot as plt

from openpyxl import load_workbook

from sklearn.preprocessing import StandardScaler

# -----------------------------------------------------------------------------
# Load the data
fnam = 'outputs_testmatrix.xlsx'
sheet = 'Test_matrix'
sheet2 = 'Univariate'

outputs = pd.read_excel(fnam,sheet,usecols='AB:BI',skiprows=[1])
sources = pd.read_excel(fnam,sheet,usecols='B:AA',skiprows=[1])
params = pd.read_excel(fnam,sheet,usecols='BJ:BQ',skiprows=[1])

# Add extra derived data to outputs
outputs['dTmax'] = outputs.apply(lambda row: row.T_C_max - row.T_G_min, axis=1)
outputs['dTavg'] = outputs.apply(lambda row: row.T_C_avg - row.T_G_avg, axis=1)
outputs['dTGmax'] = outputs.apply(lambda row: row.T_G_max - row.T_G_min, axis=1)

# Pick columns to use in the PCA
picks = [outputs.T_C_min, outputs.T_C_max, outputs.dTmax, outputs.dTavg, outputs.T_G_min, outputs.T_G_avg, outputs.dTGmax, outputs.velo_G, outputs.velo_C, outputs.K_eff, outputs.V_G, outputs.V_C ,outputs.sink_gr, sources.Qcomb, sources.m_fuel_c, sources.Qcorr, sources.Xflue, sources.S_c]
xframe = pd.concat(picks, axis=1)


# Standardise data to unit scale with mean of 0 and variance of 1
X_std = StandardScaler().fit_transform(xframe)


# Compute covariance matrix of the standardised data
xcov = np.cov(X_std.T)

# Get eigenvalues and eigenvectors for covariance matrix
eig_vals, eig_vecs = np.linalg.eig(xcov)

# Test and print to have sanity check
print('Eigenvectors \n%s' %eig_vecs)
print('\nEigenvalues \n%s' %eig_vals)

for ev in eig_vecs:
    np.testing.assert_array_almost_equal(1.0, np.linalg.norm(ev))
print('\nEverything ok!')

# Make a list of (eigenvalue, eigenvector) tuples
eig_pairs = [(np.abs(eig_vals[i]), eig_vecs[:,i]) for i in range(len(eig_vals))]

# Sort the (eigenvalue, eigenvector) tuples from high to low
eig_pairs.sort()
eig_pairs.reverse()

# Visually confirm that the list is correctly sorted by decreasing eigenvalues
print('\nEigenvalues in descending order: \n')
for i in eig_pairs:
    print(i[0])


# sum up eigenvalues to evaluate how many PC's we need.
tot = sum(eig_vals)
var_exp = [(i / tot)*100 for i in sorted(eig_vals, reverse=True)]

cum_var_exp = np.cumsum(var_exp)

# plot bar and line plot to show explained variance 

plt.figure(1,(7,5))
plt.grid(color=[0.9,0.9,0.9], linestyle='--', zorder=0)
plt.bar(np.arange(len(var_exp))+1, var_exp, label='variance explained per PC',zorder=2)
plt.plot(np.arange(len(var_exp))+1, cum_var_exp, color=[0.9,0.6,0.2], marker='.', label='cumulative explained variance')
plt.xticks(np.arange(len(var_exp))+1)
plt.yticks([0,10,20,30,40,50,60,70,80,90,100])
plt.ylabel('$\%$ explained variance')
plt.legend(loc='center right')
plt.show()





