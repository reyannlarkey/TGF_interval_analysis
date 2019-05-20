import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import socket,os
from clusterer import read_clusters
import datetime
import matplotlib
my_df = pd.read_csv("QuestTGFs_FULL_flash_reduced_No_Pulses.csv", header = None, names = ["TGF_ID", "Rey", "Will", "valid"])

TGF_list = my_df['TGF_ID'][my_df['valid']==1]

if socket.gethostname() == "reyannlarkey-personal":
    data_folder = "/home/reyannlarkey/Desktop/ENTLN_Processing/Merged/merged_clustering/"
else:
    data_folder = "/mnt/driveA/ENTLN_Processing/Merged/merged_clustering/"

def scrambler(datafolder = data_folder, scramble_n = 100000, statistic = 'Median', runs = 1, show_plot = True, save_plot = False, save_plot_thresh = 0.0):

    def calc(data, axis = 0): # Deterimines which statistic to apply to the data
        if statistic == 'median' or statistic == 'Median':
            return np.median(data, axis = axis)
        elif statistic == 'mean' or statistic == 'Mean':
            return np.mean(data, axis = axis)

    for thing in range(runs):
        time = datetime.datetime.now().time().strftime('%Hh_%Mm_%Ss')
        major_pres = np.ones((scramble_n,1))
        major_posts = np.ones((scramble_n,1))

        real_pres = []
        real_posts = []

        med_in = []
        for index, item in enumerate(TGF_list):
            print("Scrambling TGF #:", index, "Run #", thing)
            if index >=300: # This was included for testing purposes
                break
            file = data_folder+item+".txt"

            if os.path.exists(file):
                x = read_clusters(file, verbose = False)
                dt_list = x.dts


                real_pre = x.dt_pre
                real_post = x.dt_post
                med = calc(dt_list)
                real_pres.append(real_pre/med)
                real_posts.append(real_post/med)

                pres = []
                posts = []

                rand_indx = np.random.randint(low = 0, high = len(dt_list)-1, size = scramble_n)
                pre_scrambled, post_scrambled =  (dt_list[rand_indx]/med, dt_list[rand_indx+1]/med)

                pre_scrambled = pre_scrambled.reshape((scramble_n,1))
                post_scrambled = post_scrambled.reshape((scramble_n, 1))


                major_pres = np.hstack((major_pres,pre_scrambled))
                major_posts = np.hstack((major_posts, post_scrambled))
                #print(major_pres)

        major_pres = np.delete(major_pres,0,1)
        major_posts = np.delete(major_posts,0,1)

        pre_meds = calc(major_pres, axis = 1)
        post_meds = calc(major_posts, axis = 1)

        # ## COLUMN = DIFFERENT TGF!
        # ## ROW = DIFFERENT SCRAMBLE!
        # ## 206 x n_scrambles array
        #
        #
        #
        # ####################### IF YOU WANT TO PLOT ALL HISTOGRAMS! #########################
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        # bins = np.arange(0,10,0.25)
        # pre_hists = np.ones((1,len(bins)-1))
        # p2, _ = np.histogram(real_pres, bins=bins)
        #
        # for count, i in enumerate(major_pres):
        #     meds = np.median(i)
        #     plt.axvline(x= meds, color = 'black')
        #     p,_ = np.histogram(i, bins = bins)
        #     ax.step(bins[0:-1], p, where = 'post', color = 'blue')
        #
        # ax.step(bins[0:-1], p2, where = 'post', color = 'red')
        # plt.axvline(x = np.median(real_pres), color = 'red')
        # plt.xlim(0,7)
        # plt.xlabel("Pre/Median")
        # plt.show()
        # #####################################################################################
        #
        ##################### IF YOU WANT THE BOOTSTRAP RESULTS! ############################

        percent_above = (len(pre_meds[pre_meds >= calc(real_pres)])/len(pre_meds))*100.0
        percent_below = (len(post_meds[post_meds <= calc(real_posts)])/len(post_meds))*100.0
        print("Percent above {}:".format(round(calc(real_pres), 2)), percent_above)
        print("Percent below {}:".format(round(calc(real_posts), 2)), percent_below)


        pre, bins = np.histogram(pre_meds, bins = np.arange(0.5, 1.5, 0.005)+0.0025)
        post, _ = np.histogram(post_meds, bins = np.arange(0.5,1.5, 0.005)+0.0025)

        matplotlib.style.use('bmh')
        plt.figure(figsize=(13,5))
        ax = plt.subplot(111)
        real_pre_med = round(calc(real_pres),2)
        real_post_med = round(calc(real_posts),2)
        extraticks = [real_pre_med, real_post_med]
        # plt.xticks(list(plt.xticks())+extraticks)
        plt.xlim(0.7,1.3)


        lim = ax.get_xlim()
        ax.set_xticks(list(ax.get_xticks()) + extraticks)
        ax.set_xlim(lim)
        plt.step(bins[0:-1], pre, where='post', label='Pre TGF', color = '#354565', alpha = 0.75)
        plt.step(bins[0:-1], post, where='post', label='Post TGF', color = 'magenta', alpha = 0.75)
        plt.axvline(x = real_pre_med, color = '#354565', linestyle = '--')
        plt.axvline(x = real_post_med, color = 'magenta', linestyle = '--')
        plt.text(0.05, 0.5, "Percent above {}: {}, N = {}".format(round(calc(real_pres),2), round(percent_above,5), len(pre_meds[pre_meds >= calc(real_pres)])),
                 color = '#354565', transform = ax.transAxes)
        plt.text(0.05, 0.45, "Percent below {}: {}, N = {}".format(round(calc(real_posts),2), round(percent_below,5), len(post_meds[post_meds <= calc(real_posts)])),
                 color = 'magenta', transform = ax.transAxes)
        plt.legend()
        plt.title("{} Trials".format(scramble_n))
        plt.ylabel("# of Trials with Value".format(scramble_n))
        plt.xlabel("Ratio to {}".format(statistic.capitalize()))

        if show_plot == True:
            plt.show()

        if save_plot == True and percent_above > save_plot_thresh:
            plt.savefig("/home/reyannlarkey/TGF_data_analysis/Plots/BootstrappedScrambles/{}_TIME_{}_SCRAMBLES_{}_{}.pdf".format(
                                time,
                                scramble_n,
                                len(pre_meds[pre_meds >= calc(real_pres)]),
                                statistic.upper()))
        plt.close()


matplotlib.rcParams.update({'errorbar.capsize': 2})
def random_remover(percent=[0.0], ntrials = 1, statistic = 'mean'):


    def calc(data, axis=0):  # Deterimines which statistic to apply to the data
        if statistic == 'median' or statistic == 'Median':
            return np.median(data, axis=axis)
        elif statistic == 'mean' or statistic == 'Mean':
            return np.mean(data, axis=axis)


    df = pd.read_csv("232PresPostsTyp_{}.csv".format(statistic.lower()))
    pre_missing_full = {}
    post_missing_full = {}

    for percent in percent:
        pre_missing = []
        post_missing = []
        print(percent)
        for n in range(ntrials):
            remove_n =int(len(df.index) * (percent/100.0))

            drop_indices = np.random.choice(df.index, remove_n, replace=False)
            # print(drop_indices)
            df_subset = df.drop(drop_indices).reset_index(drop= True)

            pre_meds = df_subset['pres']/df_subset['{}'.format(statistic.lower())]
            post_meds  = df_subset['posts']/df_subset['{}'.format(statistic.lower())]

            pre_missing.append(calc(pre_meds))
            post_missing.append(calc(post_meds))

            # plt.scatter(df_subset['{}'.format(statistic.lower())],post_meds)

        pre_missing_full[percent] = pre_missing
        post_missing_full[percent] = post_missing

    plt.figure(figsize=(9,8))
    plt.subplot(211)
    for key, value in pre_missing_full.items():
        if np.mean(value)-3*np.std(value) <=1.0:
            col = 'red'
        else:
            col = 'gray'

        plt.plot(key, np.mean(value), marker = '.', color = '#354565')
        plt.errorbar(key, np.mean(value), yerr = 3*np.std(value),capthick=2, ecolor = col)

    plt.axhline(y = 1, color = 'black')
    plt.ylim(0,2.5)
    plt.xlim(0,100)
    plt.ylabel('Pre / {}'.format(statistic.capitalize()))
    plt.text(x = 5, y = 2, s = "Pre", color = '#354565')
    plt.text(x = 5, y = 1.75, s = "3 sigma error", color = 'gray')
    plt.title("{} Trials".format(ntrials))
    plt.subplot(212)
    for key, value in post_missing_full.items():
        if np.mean(value)+3*np.std(value) >=1.0:
            col = 'red'
        else:
            col = 'gray'
        plt.plot(key, np.mean(value), marker = '.', color = 'magenta', )
        plt.errorbar(key, np.mean(value), yerr = 3*np.std(value),capthick=2, ecolor = col)

    plt.axhline(y = 1, color = 'black')
    plt.ylim(0, 2.5)
    plt.xlim(0,100)
    plt.ylabel('Post / {}'.format(statistic.capitalize()))
    plt.xlabel('Percent of values removed from original 206')
    plt.text(x = 5, y = 2, s = "Post", color = 'magenta')
    plt.text(x = 5, y = 1.75, s = "3 sigma error", color = 'gray')
    plt.show()



# random_remover(statistic='median', percent = np.arange(0, 100, 1), ntrials = 100)




scrambler(scramble_n=5000000, save_plot=False, save_plot_thresh=0.0, show_plot=True, runs = 1000, statistic='median')