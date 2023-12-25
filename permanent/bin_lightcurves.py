import numpy as np 
import csv
import lightkurve as lk
#import exoplanet as xo
import matplotlib.pyplot as plt
import random as rd
import time
import pandas as pd
import dill 
import batman

data_points = 304128
bin_factor = 72

pd.options.display.max_columns = 10
pd.options.display.max_rows = 10

for i in range(1,25):
    camera = str(i)
    print("Computing camera", camera)
    if i < 10:
        camera = str(0)+camera
    for quarter in ["Q23", "Q24"]:
        
        with open('../jmcc/platosim/P1/per_camera_data/c'+camera+'_'+quarter+'_fluxes.pkl', 'rb') as f:
            all_lcs = dill.load(f)
        
        print("Raw lcs")
        print(all_lcs.shape)
        nb_stars = all_lcs.shape[0]
        #print(all_lcs)
        print("")
        
        binned_lcs = all_lcs.reshape(nb_stars, data_points//bin_factor, bin_factor)
        #print("Reshaped array")
        #print(all_lcs.shape)
        #print("")
        ##print(all_lcs)
        #
        binned_lcs = np.mean(binned_lcs,axis=2)
        print("Binned lcs")
        print(binned_lcs.shape)
        #print(binned_lcs)
        print("")
        print("")
        
        with open('data/binned'+str(bin_factor)+'_c'+camera+'_'+quarter+'_fluxes.pkl', 'wb') as f:
           dill.dump(binned_lcs, f)

        