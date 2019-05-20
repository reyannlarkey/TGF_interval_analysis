'''
This script will take the RAW data from a specified folder and process it.

To process data we take the following steps:

    FOR ENTLN DATA
    1. Scrape for times w/in +/- 10 mins of Fermi's reported TGF time
    2. Scrape for latitude and longitude w/in 1000km of Fermi's footprint
    3. At this point we can either merge the WWLLN Datasets, if available, or just leave the two datasets separate
    4. Cluster either the separate or merged datasets
    5. Find how many of the clusters meet the requirements we've set (i.e. the "compactness data")

I'm going to try and create a class to handle all of this
'''

import os,sys
from functools import reduce
import pandas as pd
import numpy as np
import gzip
from math import radians, cos, sin, asin, sqrt, degrees
import warnings
import clusterer
import analysis_functions as af
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore",category=DeprecationWarning)



class data_processor():

    def __init__(self, raw_data_folder = "", TGF_list =""):
        print("Running Data Processor")
        self.DataStorage = raw_data_folder
        self.TGF_list = TGF_list

    def haversine(self, lon1, lat1, lon2, lat2):
        """
        Calculate the great circle distance between two points
        on the earth (specified in decimal degrees)
        """
        # convert decimal degrees to radians
        lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

        # haversine formula
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
        c = 2 * asin(sqrt(a))
        r = 6371  # Radius of earth in kilometers. Use 3956 for miles
        return c * r



    def progressBar(self, name, value, endvalue, bar_length=70, width=30): # Because this runs pretty cool with a
                                                                           # progress bar
        percent = float(value) / endvalue
        arrow = '-' * int(round(percent * bar_length) - 1) + '>'
        spaces = ' ' * (bar_length - len(arrow))
        sys.stdout.write("\r{0: <{1}} : [{2}]{3}%".format(name, width, arrow + spaces, round(percent * 100,2)))
        sys.stdout.flush()
        if value == endvalue:
            sys.stdout.write('\n\n')

    def time_scraper(self):

        TGFList = self.TGF_list#"/home/reyannlarkey/TGF_data_analysis/UsefulTGFs_With_Dates_and_Times2.csv"
        DataStorage = self.DataStorage


        My97TGFs = pd.read_csv(TGFList, engine="python")


        # My97TGFs = My97TGFs.drop(0).reset_index()
        My97TGFs = My97TGFs.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

        d2 = pd.DataFrame(My97TGFs.UTC_Time.str.split(":").tolist(), columns="hour min sec".split())
        d2["hour"] = d2["hour"].astype(float)
        d2["min"] = d2["min"].astype(float)
        d2["sec"] = d2["sec"].astype(float)

        My97TGFs["TotalSeconds"] = d2["hour"] * 3600.0 + d2["min"] * 60.0 + d2["sec"]

        d3 = pd.DataFrame(My97TGFs.UTC_Date.str.split("-").tolist(), columns="year month day".split())

        # d3["year"] = d3["year"]
        # d3["month"] = d3["month"]
        # d3["day"] = d3["day"]

        My97TGFs['TotalDate'] = My97TGFs['UTC_Date'].str.replace("-", "")
        My97TGFs['DataFile'] = DataStorage + "LtgFlashPortions" + My97TGFs['TotalDate'] + '.csv.gz'
        My97TGFs['ENTLN_file'] = "s3cmd sync s3://eni.lightning/Archive/Pulse/" + d3['year'] + "/" + d3[
            'month'] + "/LtgFlashPortions" + My97TGFs['TotalDate'] + ".csv.gz"

        self.new_files = []
        self.new_TGF_IDs = []
        lastvalue = []
        for file in My97TGFs['DataFile']:
            if os.path.exists(file) ==True:
                p = My97TGFs['DataFile'] == file
                valid_index = My97TGFs.index[p==True]
                this_ID = My97TGFs.iloc[valid_index]['TGF_ID'].values

                if len(this_ID) == len(lastvalue) and list(this_ID) == list(lastvalue): # MAJOR PAIN IN THE BUTT
                    pass
                else:
                    lastvalue = this_ID
                    for ID in this_ID:
                        if os.path.exists(self.DataStorage+'ENTLN_Processing/TimeScraped/'+ ID + ".csv") == False: #I don't want to keep reprocessing files
                            self.new_files.append(file)
                            self.new_TGF_IDs.append(str(ID))

        self.new_TGF_IDs = np.asarray(self.new_TGF_IDs)
        print("\n Found {} Unprocessed Files \n".format(len(self.new_files)))
       # print(self.new_TGF_IDs)
        for i, filename in enumerate(self.new_files):
            self.progressBar("Time Scraping", i, len(self.new_files))
            flag = False
            TGF_ID = self.new_TGF_IDs[i]

            TGFsTotalSeconds0 = My97TGFs.query('DataFile == "{}"'.format(filename))[['TotalSeconds', 'TGF_ID']].values

            for num, time in enumerate(TGFsTotalSeconds0):
                TGFsTotalSeconds, TGF_ID = time

                flag = False
                chunksize = 10000
                count = 0

                ########################### E-N-T-L-N   F-I-L-E-S  F-O-R-M-A-T ################################
                # # ['FlashPortionID', 'FlashPortionGUID', 'FlashGUID',
                # #  'Lightning_Time', 'Lightning_Time_String', 'Latitude',
                # #  'Longitude', 'Height', 'Stroke_Type', 'Amplitude',
                # #  'Stroke_Solution', 'Offsets', 'Confidence',
                # #  'LastModifiedTime', 'LastModifiedBy']
                ###############################################################################################
                f = gzip.GzipFile(filename, 'r' )
                headerStr = S = f.readline().strip().decode().split(',')
                MyDataset = pd.DataFrame()
            #
                for chunk in pd.read_csv(filename, compression = 'gzip',  chunksize=chunksize+1, error_bad_lines=False, header= None, skiprows= 0,
                                         names = headerStr,sep=','):


                    chunk = chunk.reset_index().drop(0).reset_index()

                    chunk['Date'], chunk['Time'] = chunk['Lightning_Time_String'].str.split('T', 1).str

                    d2 = pd.DataFrame(chunk.Time.str.split(":").tolist(), columns="hour min sec".split())

                    d2["hour"] = d2["hour"].astype(float)
                    d2["min"] = d2["min"].astype(float)
                    d2["sec"] = d2["sec"].astype(float)

                    chunk["TotalSeconds"] = d2["hour"]*3600.0 + d2["min"]*60.0+ d2["sec"]

                    chunkmain = chunk[['TotalSeconds', 'Time', 'Latitude', 'Longitude', 'Height','Amplitude', 'Stroke_Type' ]]

                    p = chunkmain['TotalSeconds'].between(TGFsTotalSeconds-600.0, TGFsTotalSeconds+600.0)
                    valid_index = chunkmain.index[p==True]
                    MyDataset = pd.concat([MyDataset, chunkmain.iloc[valid_index]], ignore_index=True)

                    if len(chunkmain['TotalSeconds'][p])>=1.0:
                       # print("WE HIT PAYDIRT AT COUNT: ", count)
                        flag = True

                    if flag == True and len(chunkmain['TotalSeconds'][p])==0:
                        break #Stop seaching if we found our thing
                    if len(chunkmain['TotalSeconds']) == chunksize:
                        #print("Count = {}, TGF = {}".format(count, TGF_ID))
                        count += 1
                    chunkmain = pd.DataFrame()

                # print(MyDataset)
                MyDataset.to_csv(self.DataStorage+'ENTLN_Processing/TimeScraped/'+ TGF_ID + ".csv")

    def LatLonScraper(self, run_all = False):
        print("\n")
        TGFList_file = self.TGF_list#"/home/reyannlarkey/TGF_data_analysis/UsefulTGFs_With_Dates_and_Times2.csv"

        MyTGFs = pd.read_csv(TGFList_file, engine="python")

        if run_all == False:
            MyTGFs =MyTGFs[MyTGFs.TGF_ID.isin(self.new_TGF_IDs)].reset_index(drop =True)
        self.new_TGF_df = MyTGFs
        d2 = pd.DataFrame(MyTGFs.UTC_Time.str.split(":").tolist(), columns="hour min sec".split())
        d2["hour"] = d2["hour"].astype(float)
        d2["min"] = d2["min"].astype(float)
        d2["sec"] = d2["sec"].astype(float)

        MyTGFs["TotalSeconds"] = d2["hour"] * 3600.0 + d2["min"] * 60.0 + d2["sec"]

        km_range = 500.0

        for index, TGF in MyTGFs.iterrows():
            self.progressBar("Lat-Lon scraping", index, len(MyTGFs['TotalSeconds']))
            TGF_ID = TGF['TGF_ID'].replace("TGF", "")
            flash_file = self.DataStorage+'ENTLN_Processing/TimeScraped/'+ "{}.csv".format(TGF['TGF_ID'])
            write_file = self.DataStorage+'ENTLN_Processing/'+ 'LatLonScraped/' + TGF_ID + ".csv"
            if index < 5000 and os.path.exists(flash_file) and os.path.exists(write_file)==False:

                fermi_lat = round(float(TGF['fermi_lat']), 6)
                fermi_lon = round(float(TGF['fermi_lon']), 6)

                TGF_TIME = TGF['TotalSeconds']


                bigset = pd.DataFrame()
                chunksize = 1000
                for chunk in pd.read_csv(flash_file, chunksize=chunksize, index_col=0):
                    chunk = chunk.reset_index()
                    chunk = chunk.drop(columns=['index'])
                    valid_Idx = []
                    for i, value in enumerate(chunk['Longitude']):
                        if self.haversine(fermi_lon, fermi_lat, value, chunk['Latitude'].iloc[i]) <= km_range:
                            valid_Idx.append(i)

                    if len(valid_Idx) > 0:
                        valid = np.asarray(valid_Idx)
                        bigset = pd.concat([bigset, chunk.iloc[valid]], ignore_index=True)

                if not bigset.empty:
                    timeset = bigset['TotalSeconds'] - TGF_TIME
                    time_series = abs(timeset)
                    minindex = time_series.idxmin()

                    This_TGF = bigset.iloc[minindex].copy()

                    This_TGF['TotalSeconds'] = This_TGF['TotalSeconds'] - TGF_TIME
                    if abs(This_TGF['TotalSeconds']) <= 600.0:
                        bigset['TotalSeconds'] = round(bigset['TotalSeconds'] - TGF_TIME,6)

                        This_TGF = pd.DataFrame(data={'TotalSeconds': [bigset['TotalSeconds'].iloc[minindex]],
                                                      'Time': [bigset['Time'].iloc[minindex]],
                                                      'Latitude': [round(bigset['Latitude'].iloc[minindex],6)],
                                                      'Longitude': [round(bigset['Longitude'].iloc[minindex],6)],
                                                      'Height': [bigset['Height'].iloc[minindex]],
                                                      'Amplitude': [bigset['Amplitude'].iloc[minindex]],
                                                      'Stroke_Type': [bigset['Stroke_Type'].iloc[minindex]]})

                        DataSet = pd.concat([This_TGF, bigset], ignore_index=True)
                        DataSet.to_csv(self.DataStorage+'ENTLN_Processing/'+ 'LatLonScraped/' + TGF_ID + ".csv")

    def merger(self):
        print("\n")
        my_files = os.listdir(self.DataStorage+'ENTLN_Processing/'+ 'LatLonScraped/')
        my_files = [file.replace(".csv", "") for file in my_files]
        for i, file in enumerate(my_files):

            TGF_ID = file
            ### GET THE WWLLN DATA ###
            wl_source = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/tgf/cat_files/tgf_wwlln_data/"
            WWLLN_file = wl_source + "gbm_WWLLN.oTGF" + TGF_ID + ".txt"

            wl_data_fields = ['count', 'Longitude', 'Latitude', 'TotalSeconds', 'extras', 'extras2']
            full_WWLLN_df = pd.read_csv(WWLLN_file, skiprows=6, names=wl_data_fields, header=None, sep=r'\s*',
                                        engine='python', error_bad_lines=False)

            full_WWLLN_df['Longitude'] -= 360.00
            full_WWLLN_df['TotalSeconds'] = round(full_WWLLN_df['TotalSeconds'], 6)
            WWLLN_df = full_WWLLN_df[["Latitude", "Longitude", "TotalSeconds"]]

            ### GET THE ENTLN DATA ###
            en_source = self.DataStorage+'ENTLN_Processing/'+ 'LatLonScraped/'
            ENTLN_file = en_source + TGF_ID + ".csv"

            full_ENTLN_df = pd.read_csv(ENTLN_file, engine='python', skiprows=[1])

            # Time Correction!
            full_ENTLN_df['TotalSeconds'] = round(full_ENTLN_df['TotalSeconds'] + 0.001717, 6)
            ENTLN_df = full_ENTLN_df[['Latitude', 'Longitude', 'TotalSeconds']]
            data_frames = [WWLLN_df, ENTLN_df]
            merged_df = reduce(
                lambda left, right: pd.merge(left, right, on=['TotalSeconds', 'Latitude', 'Longitude'], how='outer'),
                data_frames)

            ### THIS IS THE MERGED DATA SET!!! ###
            merged_df = merged_df.sort_values(by=['TotalSeconds']).reset_index(drop=True)

            ### Now I want to find the TGF --> Closest absolute value to 0
            time_series = abs(merged_df['TotalSeconds'])
            minindex = time_series.idxmin()

            located_TGF = pd.DataFrame(data={'Latitude': [merged_df['Latitude'].iloc[minindex]],
                                                  'Longitude': [merged_df['Longitude'].iloc[minindex]],
                                                  'TotalSeconds': [merged_df['TotalSeconds'].iloc[minindex]]})

            final_df = pd.concat([located_TGF, merged_df], ignore_index=True)

            # writing to a folder
            write_folder = self.DataStorage+'ENTLN_Processing/'+ 'Merged/'
            final_df.to_csv(write_folder + TGF_ID + ".csv")
            self.progressBar("Merging", i, len(my_files)-1)

    def clustering(self, datakind = "MERGED"):
        print("\n")

        myTGFs= self.TGF_list#"/home/reyannlarkey/TGF_data_analysis/UsefulTGFs_With_Dates_and_Times2.csv"
        myTGFs_df = pd.read_csv(myTGFs)
        print(myTGFs_df)

        # with PdfPages('merged_maps.pdf') as pdf: # Generating a long PDF with lots of pages
        for indx, i in enumerate(myTGFs_df['TGF_ID']):
            if indx > 5000:  # This is for trouble shooting reasons.
                break
            i =i.replace("TGF", "")

            file_names = i + ".csv"

            file = self.DataStorage+'ENTLN_Processing/'+ 'LatLonScraped/' + file_names

            if os.path.exists(file):
                df = pd.read_csv(file, nrows=20)
                if len(df['TotalSeconds']) >= 5:  # I need more than 5 data points...otherwise this is dumb.

                    x = clusterer.get_clusters(tgf_ID=i, datakind = datakind)

                    ##### ONLY UNCOMMENT IF YOU WANT TO OVERWRITE YOUR CLUSTERING DATA #####
                    if datakind =="MERGED":
                        x.cluster_write(location=self.DataStorage+'ENTLN_Processing/'+ 'Merged/'+'merged_clustering/')
                    elif datakind == "ENTLN":
                        x.cluster_write(location=self.DataStorage+'ENTLN_Processing/'+'Clustered/')
                    elif datakind == "WWLLN":
                        x.cluster_write(location=self.DataStorage+'ENTLN_Processing/'+'hdbscan_clustering/')
            self.progressBar("Clustering", indx, len(myTGFs_df['TGF_ID'])-1)
                    ########################################################################


data_path = "FOLDER CONTAINING RAW DATA"
list_path = "PATH TO LIST OF TGFS WITH DATES AND TIMES"

data = data_processor(raw_data_folder=data_path, TGF_list = list_path)
print("Warning: this could take a while for a lot of files")
data.time_scraper()
data.LatLonScraper(run_all =True)
data.merger()
data.clustering(datakind = "MERGED")
# data.compactness_tester(datakind = "MERGED")