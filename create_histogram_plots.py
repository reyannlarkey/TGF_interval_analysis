'''
The purpose of this script will be to:
    (1) create the "compactness data" if it does not already exist
    (2) create the lognormal plots and histograms for the compactness data

I got tired of rewriting this all the time, so I'm going to write one script that will do it from here on out.
'''
import analysis_functions as af
import clusterer
import matplotlib.pyplot as plt
import numpy as np





#List of TGFs with Pres and posts
import pandas as pd
my_df = pd.read_csv("QuestTGFs_FULL_flash_reduced_No_Pulses.csv", header = None, names = ["TGF_ID", "Rey", "Will", "valid"])

TGF_list = my_df['TGF_ID'][my_df['valid']==1]

print("LEN = ",len(TGF_list.values))








prelist, *rest = af.logNorm_fit_and_hist(file = "QuestTGFs_FULL_flash_reduced_No_Pulses.csv", datakind = "Merged",
                                         plot = True, title = "Time Differences for: ", subs = [111], plot_KDE = False,
                                         verbose = False, statistic = "Median")


plt.show()

