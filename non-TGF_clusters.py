import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import clusterer as clusterer
from haversine import haversine
from ProgressBar import progressBar
from matplotlib.backends.backend_pdf import PdfPages

plt.style.use('bmh')


def get_max_dist(cluster_df):
    cluster_df = cluster_df[cluster_df['prob']>0.2]

    location_df = cluster_df[['lat', 'lon']]

    # print(location_df)
    dm = np.asarray([[haversine(p1[1], p1[0], p2[1], p2[0]) for p2 in location_df.values] for p1 in location_df.values])

    max_coords = np.where(dm == np.amax(dm))
    # print(np.where(dm == np.amax(dm)))


    point1 = location_df.iloc[max_coords[0][0]]
    point2 = location_df.iloc[max_coords[1][0]]

    max_dist = haversine(point1.lon, point1.lat, point2.lon, point2.lat )

    return max_dist, point1, point2



def get_valid_storms(potential_storms_file = '', valid_storms_file = ''):
    # Goes through a list of potential storms and finds only the valid ones
    df = pd.read_csv(potential_storms_file)
    MyENTLN = '/mnt/driveA/ENTLN_Processing/Merged/merged_clustering/'
    groups = df.groupby('TGF_ID')


    new_df = pd.DataFrame(columns=['TGF_ID', 'Cluster_ID'])
    for index, (name, group) in enumerate(groups):
        print(name)
        for label in group['Cluster_ID'].values:
            print(index)
            x = clusterer.read_clusters(file=MyENTLN + name + '.txt', specific_clusters=[label])

            x.get_general_flahses()
            print(set(x.clust_df['flash']))
            if max(set(x.clust_df['flash']))>13:
                max_dist, p1, p2 = get_max_dist(x.clust_df)
                if max_dist == 0:
                    clusterer.plot_indv_clusters(x, zoom=False, zoom_degree=3.5, title=x.TGF_ID, sub=111, prob_thresh=0.2, cluster_IDs=[label])
                    continue
                if max_dist<=60:
                    new_df = new_df.append({'TGF_ID':  x.TGF_ID, 'Cluster_ID': label}, ignore_index=True)
                    print(name)

            new_df.to_csv(valid_storms_file, index = False)










def scramble_flashes(valid_storms_file = "", n_storms = 0, scramble_n = 5, max_n_at_once = 10000, savefig_title = ""):
    df= pd.read_csv(valid_storms_file)
    sub_df = df.sample(n = n_storms).reset_index()



######
    dts_dict = {}
    # START BY GETTING A LIST OF ALL TIME SEPARATIONS FOR THE TGFS IN THE SUB_DF
    #----------------------------------------------------------------------------------------------
    count = 0
    groups = sub_df.groupby('TGF_ID')

    for index,  (name, group) in enumerate(groups):
        for label in group['Cluster_ID'].values: # Go get the list of dts for each of the storms!
            progressBar("Getting {} TGFS".format(n_storms), value = count, endvalue=n_storms-1)
            x = clusterer.read_clusters(file='/mnt/driveA/ENTLN_Processing/Merged/merged_clustering/'+name+'.txt', specific_clusters=[label]) #Individual Storms
            x.get_general_flahses()
            dts = x.clust_dts # time difference for all flashes in the storm
            dts_dict[count] = dts
            count+=1
    #----------------------------------------------------------------------------------------------



######
    # NOW WE NEED TO TAKE SOME NUMBER OF RANDOM VALUES FROM EACH OF THE STORMS USED!


    if max_n_at_once > scramble_n:
        max_n_at_once = scramble_n
        end_val = 1
    else:
        end_val = len(range(scramble_n//max_n_at_once))-1

    full_pre_median_list = []
    full_post_median_list = []

    for chunk in range(scramble_n//max_n_at_once):
        if end_val!=1:
            progressBar("Scrambling", value = chunk, endvalue=end_val)
        pre_stack = np.ones(max_n_at_once)
        post_stack = np.ones(max_n_at_once)

        for key, dt_list in dts_dict.items():

            random_choice_pres = np.random.choice(dt_list[0:-1],size = max_n_at_once)  # cant use the last one as a pre-interval
            random_choice_posts = np.random.choice(dt_list[1:], size = max_n_at_once)  # cant use the first one as a pre-interval

            random_choice_pres /= np.median(dt_list)
            random_choice_posts /= np.median(dt_list)


            pre_stack = np.vstack((pre_stack, random_choice_pres))
            post_stack = np.vstack((post_stack, random_choice_posts))
        #----------------------------------------------------------------------------------------------
        # Remove the 1st row because it was a dummy
        ### EACH ROW = NEW STORM
        ### EACH COLUMN = NEW RANDOM FLASH
        pre_stack = np.delete(pre_stack, 0, 0)
        post_stack = np.delete(post_stack, 0, 0)
        # Take median values
        list_of_pre_medians = np.median(pre_stack, axis=0)
        list_of_post_medians = np.median(post_stack, axis=0)
######
        full_pre_median_list = np.concatenate((full_pre_median_list,list_of_pre_medians ))
        full_post_median_list = np.concatenate((full_post_median_list, list_of_post_medians))


    if scramble_n % max_n_at_once >0:

        max_n_at_once = scramble_n % max_n_at_once
        pre_stack = np.ones(max_n_at_once)
        post_stack = np.ones(max_n_at_once)

        for key, dt_list in dts_dict.items():
            random_choice_pres = np.random.choice(dt_list[0:-1], size=max_n_at_once)  # cant use the last one as a pre-interval
            random_choice_posts = np.random.choice(dt_list[1:], size=max_n_at_once)  # cant use the first one as a pre-interval

            random_choice_pres /= np.median(dt_list)
            random_choice_posts /= np.median(dt_list)

            pre_stack = np.vstack((pre_stack, random_choice_pres))
            post_stack = np.vstack((post_stack, random_choice_posts))
        # ----------------------------------------------------------------------------------------------
        # Remove the 1st row because it was a dummy
        ### EACH ROW = NEW STORM
        ### EACH COLUMN = NEW RANDOM FLASH
        pre_stack = np.delete(pre_stack, 0, 0)
        post_stack = np.delete(post_stack, 0, 0)
        # Take median values
        list_of_pre_medians = np.median(pre_stack, axis=0)
        list_of_post_medians = np.median(post_stack, axis=0)
        print(list_of_pre_medians)
        ######
        full_pre_median_list = np.concatenate((full_pre_median_list, list_of_pre_medians))
        full_post_median_list = np.concatenate((full_post_median_list, list_of_post_medians))


    pre_hist, bins = np.histogram(full_pre_median_list, bins = np.arange(0.5, 1.5, 0.005)+0.0025)
    post_hist, _ = np.histogram(full_post_median_list, bins =np.arange(0.5, 1.5, 0.005)+0.0025)



    bin_centers =bins[0:-1]# +  0.5*(bins[1]-bins[0])

    #pdffig = PdfPages('Plots/NON_TGF_SCRAMBLES.pdf')

    plt.figure(figsize=(13, 5))
    ax = plt.subplot(111)
    plt.step(bin_centers, pre_hist,  where='post',color='#354565', label = ' \'Pres\' ')
    plt.step(bin_centers, post_hist, where='post', color = 'magenta',label = ' \'Posts\' ')

    plt.axvline(x = 1.27, color='#354565', linestyle = '--')
    plt.axvline(x = 0.97, color = 'magenta', linestyle = '--')
    plt.title( "{} Trials\nNon-TGF-Producing".format(len(full_post_median_list)))


    N_above = len(full_pre_median_list[full_pre_median_list>1.27])
    N_below = len(full_post_median_list[full_post_median_list<0.97])
    percent_above = N_above/scramble_n
    percent_below = N_below/scramble_n
    real_pre_med = 1.27
    real_post_med = 0.97

    plt.text(0.05, 0.5, "Percent above 1.27: {}, N = {}".format(round(percent_above*100,5), N_above),
             color='#354565', transform=ax.transAxes, fontsize=12)
    plt.text(0.05, 0.45,"Percent below 0.97: {}, N = {}".format(round(percent_below*100,5), N_below),
             color = 'magenta', transform=ax.transAxes,fontsize=12)

    plt.xlim(0.7, 1.3)
    extraticks = [real_pre_med, real_post_med]
    lim = ax.get_xlim()
    ax.set_xticks(list(ax.get_xticks()) + extraticks)
    ax.set_xlim(lim)
    plt.xlabel("Ratio to Median")
    plt.ylabel("# of Trials with Value")
    plt.legend(loc = 'upper left')
    if savefig_title != "":
        plt.savefig(savefig_title+'_{}.pdf'.format(N_above))
        plt.close()
    else:
        plt.show()




for i in range(1000):
    print(i, "out of 1000")
    scramble_flashes(valid_storms_file='NON_TGF_PRODUCING_SAMPLES2.csv', n_storms = 219, scramble_n= 5000000, max_n_at_once = 10000,savefig_title='Plots/Non_TGF_Scrambles/NON_TGF_PRODUCING_SCRAMBLES')

# df = pd.read_csv("NON_TGF_PRODUCING_SAMPLES2.csv")savefig_title='Plots/Non_TGF_Scrambles/NON_TGF_PRODUCING_SCRAMPLES{}.pdf'.format(i)
