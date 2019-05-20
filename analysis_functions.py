import numpy as np
from clusterer import *
from haversine import haversine  # custom made function
import cartopy.crs as ccrs
from matplotlib.gridspec import GridSpec
import pandas as pd
from matplotlib.colors import Normalize
import matplotlib.colorbar as cb
from scipy.optimize import curve_fit
import scipy.stats as stat
import os, socket
import sys

def progressBar(name, value, endvalue, bar_length=70, width=30):  # Because this runs pretty cool with a
    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length) - 1) + '>'
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write("\r{0: <{1}} : [{2}]{3}%".format(name, width, arrow + spaces, round(percent * 100, 2)))
    sys.stdout.flush()
    if value == endvalue:
        sys.stdout.write('\n\n')

from tkinter import *
from tkinter.filedialog import askopenfilename


def logNorm_fit_and_hist(file, scale="Linear", title="", prob=0.2, subs=(211, 212), orig=False,
                         plot=True, random_removal=False, removal_percent=0.0, datakind="WWLLN",
                         plot_KDE = False, verbose = False, statistic = 'median'):

    def calc(data): # Deterimines which statistic to apply to the data
        if statistic == 'median' or statistic == 'Median':
            return np.median(data)
        elif statistic == 'mean' or statistic == 'Mean':
            return np.mean(data)



    # So the user can select a file if the OS couldn't find the one they supplied
    if os.path.exists(file):
        compactness_data = pd.read_csv(file, names=['TGF_ID', 'Rey','will','TGF'], header=None)
    else:

        class MyFrame(Frame):
            def __init__(self):
                Frame.__init__(self)
                self.master.title("Oops!")
                self.master.rowconfigure(5, weight=1)
                self.master.columnconfigure(5, weight=1)
                self.grid(sticky=W + E + N + S)
                message = "Couldn't find the file, try selecting one"
                self.button = Button(self, text=message, command=self.load_file, width=len(message))
                self.button.grid(row=1, column=0, sticky=W)
                self.pathVar = StringVar()
                # self.path = Entry(textvariable=self.pathVar)
                # self.path.grid(row=1, column=1, sticky='we')

            def load_file(self):
                fname = askopenfilename(filetypes=(("CSV", "*.csv"), ("All files", "*.*")))
                if fname and os.path.exists(fname):
                    self.fname = fname
                    self.quit()

        class fileGetter():
            def __init__(self):
                file_selector = MyFrame()
                file_selector.mainloop()
                file_selector.quit()
                self.file_selected = file_selector.fname

        file = fileGetter().file_selected
        compactness_data = pd.read_csv(file, names=['TGF_ID', 'TGF'], header=None)


    items = compactness_data['TGF_ID'][compactness_data['TGF'] == 1]


    # I have three different bins here but I'm too lazy to make it cleaner
    # at some point I will go back through and clean this up.
    # Not really a priority at this point though...
    half_bins = np.arange(0.0, 600, 10.0)  # half_bins = np.logspace(0,2.778, 100)

    binwidths = np.diff(half_bins)

    # These are the lists I'm going to store data in
    datalist = np.zeros((len(half_bins) - 1))  # Histograms of original data get added to this
    postlist = []
    prelist = []
    full_dataset = []
    itemlist = []
    beforelist =[]
    afterlist = []
    avg_difflist = []
    for i, item in enumerate(items):
        if datakind == "WWLLN":
            data_folder = "hdbscan_clustering/"

        elif datakind == "ENTLN":
            if socket.gethostname() == "reyannlarkey-personal":
                data_folder = "/home/reyannlarkey/Desktop/ENTLN_Processing/Clustered/"
            else:
                data_folder = "/mnt/driveA/ENTLN_Processing/Clustered/"

        elif datakind == "MERGED" or datakind == 'Merged' or datakind =='merged':
            if socket.gethostname() == "reyannlarkey-personal":
                data_folder = "/home/reyannlarkey/Desktop/ENTLN_Processing/Merged/merged_clustering/"
            else:
                data_folder = "/mnt/driveA/ENTLN_Processing/Merged/merged_clustering/"

        file = data_folder + str(item) + '.txt'

        if os.path.exists(file):
            x = read_clusters(file, verbose = False)
            pre = x.dt_pre # Flash differences --> Not sferic difference
            post = x.dt_post # Flash differences --> Not sferic difference
            diffs = x.dts # Flash differences --> Not sferic difference


            if pre is not None and post is not None:
                # pre = abs(x.pre_flash['time_sep'].values[-1])
                # post = abs(x.post_flash['time_sep'].values[0])

                pre = x.dt_pre
                post = x.dt_post
                full_dataset.extend(diffs)
                prelist.append(pre)
                postlist.append(post)
                itemlist.append(item)
                avg_difflist.append(calc(data = diffs))

        if verbose == False:
            progressBar("Finding Pre/Post times", i, len(items)-1)

    # convert to arrays
    datalist, _ = np.histogram(full_dataset, bins = half_bins)
    full_dataset = np.asarray(full_dataset)

    prelist = np.asarray(prelist)
    postlist = np.asarray(postlist)
    datalist = np.asarray(datalist)
    avglist = np.asarray(avg_difflist)

    # medians, means, ands std. devs. of the pre and post values (of actual data, not histograms)
    med1 = calc(prelist)
    med2 = calc(postlist)
    med3 = calc(full_dataset)


    # creating the pre and post histograms
    prehist, _ = np.histogram(prelist, bins=half_bins)
    posthist, _ = np.histogram(postlist, bins=half_bins)
    avghist, _  = np.histogram(avglist, bins=half_bins)

    # print("Len of Data List", len(datalist))
    bins_to_use = half_bins[0:-1]  # Becasues np.histogram gives you 1 more value than you need
    # Taking just the values >=1.0s
    prehist_short = prehist
    posthist_short = posthist
    data_short = datalist
    binwidth = ((bins_to_use[1] - bins_to_use[0]))
    bins_short = bins_to_use#+binwidth/2

    if plot:
        ########################################################
        ### HERE'S WHERE I FIT THE DATA AND ESTIMATE ERRORS  ###
        ########################################################
        matplotlib.style.use('bmh')
        # Log Normal fitting function (Base 10 version)


        ################################
        ### PLOTTING HISTOGRAM DATA  ###
        ################################
        if scale == 'Log' or scale == "log":
            text_xpos = 0.05
            text_ypos = 0.7
        else:
            text_xpos = 0.5
            text_ypos = 0.5


        left, width = 0.05, 0.3
        bottom, height = 0.1, 0.75
        left_h = left + width + 0.01
        plt.figure(figsize=(15, 7))
        #plt.suptitle("Pre TGF\n{} Data".format(datakind), fontsize = 20)
        rect_scatter = [left, bottom, width, height]

        rect_histy = [left_h, bottom, 0.12, height]

        axScatter = plt.axes(rect_scatter)
        axScatter.axhline(y = 1, color = 'black')
        axScatter.scatter(avg_difflist, prelist / avg_difflist, color='#354565')#, label='Pre-TGF/{}'.format(statistic.capitalize()))
        axScatter.axhline(y=calc(prelist / avg_difflist), color='cyan', label='{} of Ratios = {}'.format(statistic.capitalize(), round(calc(prelist/avg_difflist),2)))
        axScatter.legend(prop={'size': 18},loc = 'upper left')
        axScatter.set_xlim(0)
        axScatter.set_ylim(0, 7)
        axScatter.set_ylabel("Pre-TGF Interval / {} Interval".format(statistic.capitalize()),fontsize = 20)
        axScatter.set_xlabel("{} inter-flash interval (s)".format(statistic.capitalize()), fontsize = 20)
        axHisty = plt.axes(rect_histy)
        axHisty.axhline(y = 1, color = 'black')
        axHisty.hist(prelist/avg_difflist, bins =np.arange(0,10,0.25), orientation = 'horizontal', facecolor = '#354565')
        axHisty.set_ylim(0,7)
        axHisty.axhline(y=calc(prelist / avg_difflist), color='cyan', label = "{} = {}".format(statistic.capitalize(), round(calc(prelist / avg_difflist), 2)))
        #axHisty.legend(prop={'size': 12}, loc = 'upper left')
        from matplotlib.ticker import NullFormatter
        nullfmt = NullFormatter()
        axHisty.yaxis.set_major_formatter(nullfmt)

        left2, width2 = 0.53, 0.3
        left_h2 = left2 + width2 + 0.01
        # plt.figure(figsize=(9, 8))
        plt.suptitle("Pre and Post-TGF Intervals\n{} Data".format(datakind), fontsize=20)
        rect_scatter = [left2, bottom, width2, height]
        rect_histy = [left_h2, bottom, 0.12, height]

        axScatter = plt.axes(rect_scatter)
        axScatter.axhline(y=1, color='black')
        axScatter.scatter(avg_difflist, postlist / avg_difflist, color='magenta')#, label='Post-TGF/{}'.format(statistic.capitalize()))
        axScatter.axhline(y = calc(postlist/avg_difflist), color = 'cyan', label = '{} of Ratios = {}'.format(statistic.capitalize(), round(calc(postlist/avg_difflist),2)))
        axScatter.legend(prop={'size': 18}, loc = 'upper right')
        axScatter.set_xlim(0)
        axScatter.set_ylim(0, 7)
        axScatter.set_ylabel("Post-TGF Interval / {} Interval".format(statistic.capitalize()), fontsize = 20)
        axScatter.set_xlabel("{} inter-flash interval (s)".format(statistic.capitalize()), fontsize = 20)
        axHisty = plt.axes(rect_histy)
        axHisty.axhline(y=1, color='black')
        axHisty.axhline(y=calc(postlist / avg_difflist), color='cyan', label = "{} = {}".format(statistic.capitalize(), round(calc(postlist/avg_difflist),2)) )
        axHisty.hist(postlist / avg_difflist, bins=np.arange(0, 10, 0.25), orientation='horizontal', facecolor = 'magenta' )
        axHisty.set_ylim(0, 7)
        #axHisty.legend(prop={'size': 12}, loc = 'upper left')
        from matplotlib.ticker import NullFormatter
        nullfmt = NullFormatter()
        axHisty.yaxis.set_major_formatter(nullfmt)


        
        if plot_KDE:
            from sklearn.neighbors.kde import KernelDensity
            print(avg_difflist)
            list2 = [prehist_short, posthist_short, avghist]#, data_short]
            list3 = [prelist, postlist, avg_difflist]#,full_dataset]


            colorlist = ['#345565', 'magenta','red']
            labellist =['pre','post', 'full']

            colorlist = ['blue', 'magenta','red']
            labellist =['pre','post', 'typical (MEAN)']

            sublist = [411,412,413, 414]
            fig = plt.figure(figsize=(9,8))
            fig.suptitle(title + " " + str(len(prelist)) + " TGF events")
            #plt.tight_layout()
            for indx, data in enumerate([prelist, postlist, full_dataset, []]):
                if sublist[indx] == 414:
                    plt.subplot(414)
                    for indx2, thing in enumerate(list2):
                        # kde = KernelDensity(kernel='tophat', bandwidth=10).fit(np.asarray(thing).reshape(-1, 1))
                        # s = np.linspace(0, 600, 5000)
                        # e = kde.score_samples(s.reshape(-1, 1))
                        plt.plot(bins_short+binwidth/2, thing / sum(thing), color=colorlist[indx2], linestyle='-', label = labellist[indx2])

                        plt.legend()
                        plt.xlabel("Time Sep. (s)")
                        plt.ylabel("Density")
                        plt.xlim(0,310)
                else:
                    if sublist[indx] == 411:
                        plt.title(title + " " + str(len(prelist)) + " TGF events")

                    ax = plt.subplot(sublist[indx])


                    plt.bar(bins_short+binwidth/2, list2[indx], width = binwidth, facecolor = colorlist[indx], label=labellist[indx])
                    plt.text(text_xpos, text_ypos, r"$med = $ {}".format(round(np.median(list3[indx]), 2)), transform=ax.transAxes, color=colorlist[indx])
                
                    plt.legend()
                    plt.ylabel("Events/10s")
                    plt.xlim(0,310)
            plt.show()
        return (0, 0, 0, 0, 0,0)
    else:
        print("Not plotting")
        prelist_short = prelist[prelist >= 1.0]
        postlist_short = postlist[postlist >= 1.0]
        datalist_short = full_dataset


        # print("Mean of Prelist  = {}".format(mean1))
        # print("Mean of Postlist = {}".format(mean2))
        print("Median of Prelist = {}".format(med1))
        print("Median of Postlist = {}".format(med2))
        return (prelist_short, postlist_short, datalist_short,beforelist, afterlist)



