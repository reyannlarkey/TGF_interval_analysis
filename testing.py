#import clusterer
import analysis_functions as af
import os
import pandas as pd
#import side_by_side_maps
import matplotlib.pyplot as plt
import numpy as np
import socket
import hdbscan
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as mpatches


#
if __name__ == "__main__":
    ### For the full dataset
    if socket.gethostname() == "reyannlarkey-personal":
        MyENTLN = "/home/reyannlarkey/Desktop/ENTLN_Processing/Merged/merged_clustering/"

    else:
        MyENTLN = '/mnt/driveA/ENTLN_Processing/Merged/merged_clustering/'



    ### For doing the narrowed down version
    # df = pd.read_csv("QuestTGFs_FULL.csv", header = None, names = ["TGF_ID", "valid"])
    # scraped_df = df['TGF_ID'][df['valid']==1]
    # files = MyENTLN+scraped_df+".txt"
    # files = [x for x in files.values if os.path.exists(x)]
    #
    ### For Testing
    test_dir = MyENTLN
    test_file = [MyENTLN + '100304844.txt']
    test_list = [MyENTLN + x for x in ['120814947.txt', '120706651.txt', '101008251.txt']]
    files_in_test_dir = [MyENTLN+x for x in sorted(os.listdir(test_dir))]
#
    for files in test_list:
        x = clusterer.read_clusters(file = files)
        plt.figure()
        clusterer.plot_indv_clusters(x, zoom = True, title = "TGF " + x.TGF_ID)

        # clusterer.plot_indv_clusters(x, zoom = False, title="TGF " + x.TGF_ID)
        # plt.show()
    plt.show()

    #

    # df = pd.read_csv("QuestTGFs_FULL_flash_reduced_With_CLASS.csv", header = None, names = ["TGF_ID", "Rey", "Will", "Valid", "Class"])
    # counts = (df['Class'].value_counts())
    #
    # print(counts/counts.sum() *100.0)







#
#     prelist = []
#     postlist = []
#     for file in files:
#         print(file)
#         cluster = clusterer.read_clusters(file = file, prob_thresh=0.2)
#         if cluster.pre_flash is not None:
#             prelist.append(cluster.TGF_time - (cluster.pre_flash['time_sep']))
#             postlist.append((cluster.post_flash['time_sep'])-cluster.TGF_time)
#         # print(prelist)
#     pre,bins = np.histogram(prelist, bins = np.arange(1.0,600,10))
#     posts, _ = np.histogram(postlist, bins = np.arange(1.0,600,10))
#     binwidth = (bins[1]-bins[0])/2
#     bins = bins[0:-1]+binwidth
#
#     plt.figure()
#
#     plt.bar(bins, posts, width=bins[1] - bins[0], facecolor="None", edgecolor = "magenta", picker = True)
#     plt.bar(bins, pre, width = bins[1]-bins[0], facecolor = "None", edgecolor = "blue", picker = True)
#
#     # plt.plot(prelist, postlist, 'b.')
#     plt.show()
