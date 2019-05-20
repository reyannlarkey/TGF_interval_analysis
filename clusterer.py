import numpy as np
import pandas as pd
import hdbscan
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib.colorbar as cb
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import itertools
# import seaborn as sns


import matplotlib
#matplotlib.style.use('dark_background')


useful_tgfs = "/home/reyannlarkey/Desktop/TGF_Data_Analysis/flash_rate/useful_tgfs.txt"
source = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/tgf/cat_files/tgf_wwlln_data/"


#Actually performs the clustering algorithms and saves that data to files
class get_clusters:
    def __init__(self, tgf_ID = "", datakind = "WWLLN"):

        if tgf_ID == "":
            print("NO TGF")
        else:
            self.TGF_ID = tgf_ID
            if datakind == "WWLLN":
                url = source + "gbm_WWLLN.oTGF" + tgf_ID +".txt"
                tgf_fields = ['num','Longitude', 'Latitude', 'TotalSeconds','SIMUs', 'extras']
                self.TGF_number = pd.read_csv(url, skiprows=4, sep=r'\s*', nrows=1, header=None, engine='python',
                                              names=tgf_fields)['num'].values
                self.TGF_info = tgf_info = pd.read_csv(url, skiprows=5, names=tgf_fields, header=None, sep=r'\s*',
                                                       engine='python', nrows=1)
                data_fields = ['count', 'Longitude', 'Latitude', 'TotalSeconds','extras', 'extras2']
                df = pd.read_csv(url, skiprows = 6, names = data_fields, header = None, sep=r'\s*',
                             engine ='python',  error_bad_lines=False)

                self.TGF_lat = self.TGF_info['Latitude'].values
                self.TGF_lon = self.TGF_info['Longitude'].values
            if datakind == "ENTLN":
                url = "/mnt/driveA/ENTLN_Processing/LatLonScraped/"+tgf_ID+".csv"
                self.TGF_info = pd.read_csv(url, nrows = 1, index_col=0)
                df = pd.read_csv(url, skiprows=[0,1], engine='python', error_bad_lines=False, header = None,
                                 names = ['TotalSeconds', 'UTC','Latitude', 'Longitude', 'alt','SIMUs', 'extras'])
                df =df.reset_index()
                self.TGF_lat = self.TGF_info['Latitude'].values[0]
                self.TGF_lon = self.TGF_info['Longitude'].values[0]

            if datakind == "MERGED":
                url = "/mnt/driveA/ENTLN_Processing/Merged/" + tgf_ID + ".csv"
                self.TGF_info = pd.read_csv(url, nrows=1, index_col=0)
                df = pd.read_csv(url, skiprows=[0, 1], engine='python', error_bad_lines=False, header=None,
                                 names=['Latitude', 'Longitude','TotalSeconds'])
                df = df.reset_index()
                self.TGF_lat = self.TGF_info['Latitude'].values[0]
                self.TGF_lon = self.TGF_info['Longitude'].values[0]


            self.TGF_time = self.TGF_info['TotalSeconds'].values[0]

            # Get all the other Data


            location_df = pd.DataFrame({
                'lat': df['Latitude'].values,
                'lon': df['Longitude'].values
            })


            self.location_df = pd.DataFrame({
                'lat': df['Latitude'].values,
                'lon': df['Longitude'].values
            })


            self.time_df = pd.DataFrame({
                't' : df['TotalSeconds'].values
            })

            # cluster the data
            clusterer = hdbscan.HDBSCAN(metric = 'haversine',cluster_selection_method='eom',
                                        allow_single_cluster=True).fit(np.radians(location_df))
            #plt.figure()
            #clusterer.condensed_tree_.plot(select_clusters=True, selection_palette=sns.color_palette())
#
            self.cluster_probs = clusterer.probabilities_
            self.cluster_labels = clusterer.labels_
            self.outlier_scores = clusterer.outlier_scores_

            # pick out which cluster the TGF is into
            self.TGF_index = (df[df['TotalSeconds']==self.TGF_time].index).values

            self.TGF_cluster  = self.cluster_labels[self.TGF_index][0]

            #
            self.TGF_prob = self.cluster_probs[self.TGF_index][0]
            self.TGF_outlier = self.outlier_scores[self.TGF_index][0]


    def cluster_write(self, location = ""):
        if location == "":
            print("Cant write file: No file location specified")

        else:
            #response = input("You are about to write some data, are you sure you want to proceed? (y/n)")

            #if response == "Y" or response == "y":
            f = location + self.TGF_ID +'.txt'
            with open(f, 'w+') as file:
                file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(self.TGF_lat, self.TGF_lon,
                                                             self.TGF_time, self.TGF_cluster, self.TGF_prob,
                                                             self.TGF_outlier))

                for index, line in enumerate(self.location_df['lat'].values):
                    file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(self.location_df['lat'].ix[index],
                                                                 self.location_df['lon'].ix[index],
                                                                 self.time_df['t'].ix[index], self.cluster_labels[index],
                                                                 self.cluster_probs[index], self.outlier_scores[index]))

            #print(self.TGF_ID+ ".txt has been written")

            # else:
            #     print("Not writing the data file")



# Reads in the clustering data from the text files that "get_clusters" creates
class read_clusters:
    def __init__(self, file = "", prob_thresh =0.2, verbose = False):
        self.verbose = verbose
        if file == "":
            if self.verbose:
                print("Can't read file, no file specified")
        else:
            splitfile = file.split('/')
            self.TGF_ID = splitfile[-1].split('.')[0]
            fields = ['lat', 'lon', 'time_sep', 'cluster','prob','outlier']
            tgf_line = pd.read_csv(file, nrows=1, sep='\t', header=None, names=fields, engine = 'python')

            self.TGF_line = tgf_line
            self.TGF_lat = np.round(tgf_line['lat'].values,6)

            self.TGF_lon = np.round(tgf_line['lon'].values,6)
            self.TGF_time = np.round(tgf_line['time_sep'].values,6)
            self.TGF_cluster  = tgf_line['cluster'].values
            self.TGF_prob = tgf_line['prob'].values
            self.TGF_outlier = tgf_line['outlier'].values

            df = pd.read_csv(file, sep='\t', header=None, skiprows=1, names=fields)
            self.Npts = len(df['lat'])
            self.df = df
            self.location_df = pd.DataFrame({
                'lat': df['lat'].values,
                'lon': df['lon'].values
            })

            self.cluster_probs = df['prob'].values
            self.cluster_labels = df['cluster'].values
            self.lats = df['lat'].values
            self.lons = df['lon'].values
            self.times  = df['time_sep'].values
            self.outliers = df['outlier'].values

            rdata = {}
            rdata['lats'] = {}
            rdata['lons'] = {}
            rdata['times'] = {}
            rdata['probs'] = {}

            for i, item in enumerate(self.cluster_labels):
                mask = np.nonzero(np.isin(self.cluster_labels, item))
                rdata['lats'][item] = np.take(self.lats, mask)[0]
                rdata['lons'][item] = np.take(self.lons, mask)[0]
                rdata['times'][item] = np.take(self.times, mask)[0]
                rdata['probs'][item] = np.take(self.cluster_probs, mask)[0]
            self.rdata = rdata
            self.cluster_IDs = np.asarray(list(set(self.cluster_labels)))

            self.get_pre_post_times(prob_thresh=prob_thresh)
            self.get_flashes(prob_thresh=prob_thresh)

    def get_pre_post_times(self, prob_thresh = 0.2):
        TGF_df = self.df[self.df['cluster'] == self.TGF_cluster[0]]
        TGF_df = TGF_df[TGF_df['prob'] > prob_thresh]
        # print(TGF_df)

        self.TGF_cluster_len  = len(TGF_df['lat'])
        location = TGF_df.loc[(round(TGF_df['time_sep'],6) == self.TGF_time[0]) & round(TGF_df['lat'],6).isin(self.TGF_lat)].index.values
        if len(location) == 0:
            if self.verbose:
                print("TGF sferic not associated with cluster at this probability threshold", self.TGF_ID)
            self.pre_flash, self.post_flash = (None, None)
            self.TGF_flash= 0
            self.dts =  None
            self.dt_pre = None
            self.dt_post= None
            self.TGF_flash_full = None

        else:
            TGF_index = location[0]

            pre_TGF_df = TGF_df.loc[:TGF_index].reset_index(drop = True)
            post_TGF_df = TGF_df.loc[TGF_index:].reset_index(drop = True)

            if len(post_TGF_df[post_TGF_df['time_sep']>1.0]['time_sep'])>=1 and len(pre_TGF_df[pre_TGF_df['time_sep']< -1.0]['time_sep'])>=1:
                self.post_flash = (post_TGF_df[post_TGF_df['time_sep']>1.0].iloc[0])
                self.pre_flash = (pre_TGF_df[pre_TGF_df['time_sep'] < -1.0].iloc[-1])
            else:
                if self.verbose:
                    print("Couldn't find a pre/post flash for this cluster", self.TGF_ID)
                self.pre_flash, self.post_flash = (None, None)
                self.TGF_flash = 0
                self.dts = None
                self.dt_pre = None
                self.dt_post = None
                self.TGF_flash_full = None

    def get_flashes(self, prob_thresh = 0.2):
        TGF_df = self.df[self.df['cluster'] == self.TGF_cluster[0]]
        TGF_df = TGF_df[TGF_df['prob'] > prob_thresh]
        #print(TGF_df)
        if len(TGF_df)>1 and self.pre_flash is not None and self.post_flash is not None:
            # I want to group things by times
            from sklearn.neighbors.kde import KernelDensity
            kde = KernelDensity(kernel='gaussian', bandwidth=0.25).fit(np.asarray(TGF_df['time_sep']).reshape(-1,1))
            s = np.linspace(-1000,1000,5000)
            e = kde.score_samples(s.reshape(-1,1))
            # plt.plot(s, np.exp(e))
            # plt.xlim(-600,600)
            # plt.ylim(0,0.05)
            #
            # plt.plot(TGF_df['time_sep'],np.zeros(len(TGF_df['time_sep'])),'b*')
            #


            from scipy.signal import argrelextrema
            mi, ma = argrelextrema(e, np.less)[0], argrelextrema(e, np.greater)[0]
            cuts = s[mi]

            cuts = np.insert(cuts, 0, -1000.0)
            cuts = np.append(cuts, 1000.0)

            flash_list = []
            #print(self.TGF_ID)
            for indx, cut in enumerate(cuts):
                if indx == len(cuts)-1:
                    break
                group = TGF_df[TGF_df['time_sep'].between(cut, cuts[indx+1])]
                for _ in group['time_sep']:
                    flash_list.append(indx)

            TGF_df['flash'] = flash_list
            location = TGF_df.loc[(round(TGF_df['time_sep'],6) == round(self.TGF_time[0],6)) & round(TGF_df['lat'],6).isin(self.TGF_lat)].index.values

            self.TGF_flash = TGF_df.ix[location]['flash'].values[0]
            self.TGF_flash_full = TGF_df[TGF_df['flash'] == self.TGF_flash]
            self.TGF_df = TGF_df
            self.pre_flash = self.TGF_df[self.TGF_df['flash']==self.TGF_flash -1]
            self.post_flash = self.TGF_df[self.TGF_df['flash'] == self.TGF_flash +1]
            if len(self.post_flash['flash']) ==0:
                self.pre_flash, self.post_flash = (None, None)
                self.TGF_flash = 0
                self.dts = None
                self.dt_pre = None
                self.dt_post = None
                self.TGF_flash_full = None
            #plt.plot(cuts,np.zeros(len(cuts)), 'r.')
            #plt.show()
            else:
                dts = []
                for flash in set(TGF_df['flash']):
                    if flash > 0:
                        this_flash_chunk = TGF_df[TGF_df['flash'] == flash]
                        this_flash_start = this_flash_chunk.iloc[0]['time_sep']

                        prev_flash_chunk = TGF_df[TGF_df['flash'] == flash - 1]
                        prev_flash_end = prev_flash_chunk.iloc[-1]['time_sep']

                        dt = this_flash_start - prev_flash_end
                        dts.append(dt)
                self.dts = np.asarray(dts)
                self.dt_pre = self.TGF_flash_full['time_sep'].values[0] - self.pre_flash['time_sep'].values[-1]
                self.dt_post = self.post_flash['time_sep'].values[0] - self.TGF_flash_full['time_sep'].values[-1]
        else:
            TGF_df['flash'] = 0
            self.TGF_df = TGF_df



def plot_indv_clusters(cluster, cluster_IDs = [], prob_thresh = 0.2, zoom = False, zoom_degree = 1.0, title = "", sub = 111,
                       optional_inds  = [], optional_text = "", optional_text_loc = [], TGF_legend_loc = [0.05, 0.9],
                       TGF_star_color = "lime", opt_text_color = "blue"):
    if cluster_IDs == [] or cluster_IDs== 'all':
        cluster_IDs = cluster.cluster_IDs
    #plt.subplot(sub)
    rdata = cluster.rdata

    color_map = cm.jet(np.linspace(0.0, 1.0, max(set(cluster.cluster_labels)) + 2))
    ax = plt.subplot(sub, projection=ccrs.PlateCarree())
    ax.coastlines('50m')
    land_50m = cfeature.NaturalEarthFeature('physical','land', '50m',
                                            facecolor=cfeature.COLORS['land'])
    ax.add_feature(land_50m)
    ax.text(-0.125, 0.55, 'Latitude', va='bottom', ha='center',
            rotation='vertical', rotation_mode='anchor',
            transform=ax.transAxes)
    ax.text(0.5, -0.1, 'Longitude', va='bottom', ha='center',
            rotation='horizontal', rotation_mode='anchor',
            transform=ax.transAxes)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = True
    gl.ylabels_right = False
    plt.title(title)

    props = dict(boxstyle='round', facecolor='white', alpha=0.7)
    opt_props = dict(boxstyle='round', facecolor='white', alpha=0.7)
    if len(optional_inds)>1:
        print(cluster.location_df.iloc[optional_inds[0]]['lon'])
        print(cluster.location_df.iloc[optional_inds[0]]['lat'])
        ax.scatter(cluster.location_df.iloc[optional_inds[0]]['lon'], cluster.location_df.iloc[optional_inds[0]]['lat'],
                   color='blue', edgecolor='black', marker='*', s = 100, zorder = 1000)
        ax.scatter(cluster.location_df.iloc[optional_inds[1]]['lon'], cluster.location_df.iloc[optional_inds[1]]['lat'],
                   color='magenta', edgecolor='green', marker='*', s = 100,zorder=1000)
    if optional_text != "":
        if optional_text_loc != []:
            plt.text(optional_text_loc[0], optional_text_loc[1], optional_text,
                     horizontalalignment='left',
                     verticalalignment='center', fontsize=10,
                     transform=ax.transAxes,
                     color=opt_text_color, bbox=opt_props)
        else:
            print("No specified text location, not adding text")

    marker = itertools.cycle(('.', 'v', '+', 'd', '1', 'x'))
    for cluster_ID in cluster_IDs:
        if cluster_ID == -1:
            marker = itertools.cycle(('.'))

        if cluster_ID == cluster.TGF_cluster:
            if TGF_legend_loc is not None:
                plt.text(TGF_legend_loc[0],TGF_legend_loc[1], 'TGF Lat: {} \nTGF Lon: {}'.format(cluster.TGF_lat[0], cluster.TGF_lon[0]),
                         horizontalalignment='left',
                         verticalalignment='center', fontsize=10,
                         transform=ax.transAxes,
                         color='black', bbox=props)
            #
            if cluster.TGF_prob >= prob_thresh:
                ax.scatter(cluster.TGF_lon, cluster.TGF_lat, color = TGF_star_color, edgecolor='black', marker='*', s=100, alpha=1,
                          zorder=1000)
            else:
                ax.scatter(cluster.TGF_lon, cluster.TGF_lat, color='black', edgecolor='black', marker='*', s=100,
                           alpha=1,
                           zorder=1000)
        filt = rdata['probs'][cluster_ID] >= prob_thresh
        filt2 = rdata['probs'][cluster_ID] < prob_thresh
        if len(rdata['lons'][cluster_ID][filt])>0:
            mark = next(marker)
            
            ax.scatter(rdata['lons'][cluster_ID][filt], rdata['lats'][cluster_ID][filt],
                       color=color_map[cluster_ID], marker=mark, s=20, zorder = 100)
            ax.plot(rdata['lons'][cluster_ID][filt2], rdata['lats'][cluster_ID][filt2],
                       markerfacecolor="None", markeredgecolor = 'black', markeredgewidth = 1, marker=mark, linestyle = "None", markersize =5)

        elif len(rdata['lons'][cluster_ID][filt2]) >0:
            ax.plot(rdata['lons'][cluster_ID][filt2], rdata['lats'][cluster_ID][filt2],
                    markerfacecolor="None", markeredgecolor='black', markeredgewidth=1, marker=next(marker),
                    linestyle="None",zorder = 10, markersize = 5)

        else:
            pass

    if zoom:
        ax.set_extent([cluster.TGF_lon - zoom_degree, cluster.TGF_lon +zoom_degree, cluster.TGF_lat - zoom_degree, cluster.TGF_lat + zoom_degree])
    else:
        ax.set_extent([-180, 180, -90, 90])
    return ax

