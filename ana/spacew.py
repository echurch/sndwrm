from ROOT import TFile, gPad, gStyle, TH1F, TH2F, TGraph
import pandas as pd
import pdb
import glob
import numpy as np
from matplotlib import pyplot as plt


#file_path="../data/job6_EE.xlsx"
#file_path="../data/job6b.xlsx" # use hi
file_path="../data/job1c.xlsx" # use li

df_li = pd.read_excel(file_path, sheet_name="li")
#df_hi = pd.read_excel(file_path, sheet_name="hi")
dfs = [df_li]
#dfs = [df_hi]
root_file = "./tgraph-spectra_cosmic/job6b_spectra.root"
root_file = "./tgraph-spectra_solar/job1c_spectra.root"
f_p = TFile(root_file, "RECREATE")

for df in dfs:

    # Assuming the first row is the common X-axis values
    x_values = np.array(df.columns)[1:]
    x_values = x_values.astype(float)

    xcenters = x_values
    all_data = df.iloc[1:, :]

    npts = len(x_values)

    for jj in range(all_data.shape[0]): # loop through rows
        gname = all_data.iloc[jj,0]
        title = "spectra"

        gpart = TGraph(npts,np.array(xcenters),np.array(all_data.iloc[jj,1:].values.astype(float)))
        gpart.SetTitle(gname+"_"+title)
        gpart.SetName(gname+"_"+title)
        gpart.SetMarkerStyle(20+jj)
        gpart.Write()

        filename = "./tgraph-spectra_solar/"+gname+"_"+title+".txt"
        with open(filename, 'w') as fname:
            for ii in range(npts):
                gpart.GetPoint(ii,xcenters[ii],all_data.iloc[jj,1:].values.astype(float));
                fname.write(f"{xcenters[ii]}  {all_data.iloc[jj,1:].values.astype(float)[ii]}\n")

        del gpart

    
f_p.Close()
