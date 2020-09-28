#!/usr/bin/env python
# coding: utf-8

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import os
from pathlib import Path
from math import *
import datetime
from slope_analysis import *

#########!!!INPUT!!!########!!!INPUT!!!####!!!INPUT!!!##################
#########!!!INPUT!!!########!!!INPUT!!!####!!!INPUT!!!##################
#Sample ID

#edited by Simon Peterson for Katey Anthony and Nick Hasson
#08/20/2020

#TODO- make sample ID print to command line, automatically create the ID tag
#string names of the valid measurement devices:
measurement_devices = ["bucket_sediment","bucket_snow", "chamber"]
#greenhouse gases that can be measured:
gases = ["CO2","CH4"]
#current working directory
cwd = os.getcwd()
output_folder = os.path.join(cwd,"outputs")
lgr_input_folder = os.path.join(cwd, "inputs","lgr")
weather_input_folder = os.path.join(cwd, "inputs", "weather_data") 
print("pulling lgr data from: \n" + lgr_input_folder)
print("pulling weather data, when not from internet, from: \n" + weather_input_folder + "\n")

#the location of the r_2 excel file
r_2_values_path = os.path.join(cwd, "data","r_values.xlsx")
r_2_values = pd.read_excel(r_2_values_path, header = 0)
r_2_values = r_2_values['r_values'].values.tolist()
print(r_2_values)

#read the excel file for the input. Close the file before reading it.
master_csv_path = os.path.join(cwd, "data","simon_masters.xlsx")
master_data = pd.read_excel(master_csv_path, header = 0) #headers are the first row
#apparently excel formats their data with the dates also including a timestamp. Annoying
master_data["date_(yyyy-mm-dd)"] = master_data["date_(yyyy-mm-dd)"].apply(lambda x: x.date().isoformat())


print(master_data["date_(yyyy-mm-dd)"])

#temperature and pressure data
weather_input_filepath = os.path.join(weather_input_folder, 'June2019_T_P.csv')
t_p_data = pd.read_csv(weather_input_filepath, delimiter=',', parse_dates=[['date','time']])
t_p_data['date_time'] = pd.to_datetime(t_p_data['date_time'])
#P_Pa = dfp.loc[dfp['date_time'] == '6/11/2019 0:00','air_p_mean_Pa'].values
print("using as master spreadsheet:\n" + master_csv_path)
print("The master data sheet:")
print(master_data)


#read the data from the lgr files. Clean any "dirty" files.
#TODO- make so that junk at the end of the files is cleaned up
first = True
for file in os.listdir(lgr_input_folder):
	print(file)
	if file.endswith(".txt"):
		with open(os.path.join(lgr_input_folder,file),'r') as f:
			file_text = f.read()
			if "BEGIN PGP MESSAGE" in file_text:
				print("found text \" BEGIN PGP MESSAGE \", will clean file")
				f.close()
				with open(os.path.join(lgr_input_folder, file),'r') as f:
					lines = f.readlines()
					f_new = open(os.path.join(lgr_input_folder, file.strip(".txt")) + "_cleaned.txt",'w+')
					for line in lines:
						if "-----BEGIN PGP MESSAGE-----" not in line:
							f_new.write(line)
						else:
							break
					print("successfully cleaned " + file)
					f.close()
					os.remove(os.path.join(lgr_input_folder,file))
					file = file.strip(".txt") + "_cleaned.txt"
					print("removed old file, created " + file)
					f_new.close()
		if first:
			lgr_data = pd.read_csv(os.path.join(lgr_input_folder,file), delimiter=',', header = 1, index_col = 1)
			first = False
		else:
			print("appending data")
			new_data = pd.read_csv(os.path.join(lgr_input_folder,file), delimiter=',', header = 1, index_col = 1)
			lgr_data = lgr_data.append(new_data)
		print("size of total LGR data array: " + str(lgr_data.shape))
			
#TODO- print the total lgr data into a nice csv :)
#Make the lgr date index with date and time
lgr_data.index = pd.DatetimeIndex(lgr_data.index)


#now we will create a new array with the times that we would like to add to the list
#make a for loop to run through the data which has not already been run by the program
torun_rows = []
for row in range(master_data.shape[0]):
	print(master_data.iloc[row]['program_run?'])
	if master_data.iloc[row]['program_run?'] != 'y':
		torun_rows.append(row)
print("the rows that will be run, starting at index 0:")
print(torun_rows)
print("will run the program for a total of " + str(len(torun_rows)) + " times")


#dictionary with the rows as keys and the sample IDs as values
row_ID = {}
row_gases = {}
for row in torun_rows:
	#create the sample ID. The format is:
	#yyyy-mm-dd_hh:mm:ss_location_collection-instrument
	#where the hh:mm:ss is the START TIME
	sample_ID = ""
	if master_data.iloc[row]["date_(yyyy-mm-dd)"] != np.nan:
		sample_ID =  str(master_data.iloc[row]["date_(yyyy-mm-dd)"]).replace('-','_')
	else:
		print("need a date for row: " + str(row))
		exit()
	if isinstance(master_data.iloc[row]["start_time_(hh:mm:ss)"], datetime.time):
		sample_ID = sample_ID + "_" + str(master_data.iloc[row]["start_time_(hh:mm:ss)"]).replace(":","h",1).replace(":",'m',1) + "s"
	else:
		print("no start time found for row: " + str(row))
		exit()
	if not isinstance(master_data.iloc[row]["stop_time_(hh:mm:ss)"], datetime.time):
		print("no stop time found for row: " + str(row))
		exit()
	if master_data.iloc[row]["location_(lake)"] != np.nan:
		sample_ID = sample_ID + "_" + str(master_data.iloc[row]["location_(lake)"])
	else:
		print("no location for row: " + str(row))
		sample_ID = sample_ID + "_unentered"
	if master_data.iloc[row]["measurement_device"] not in measurement_devices:
		print("invalid measurement_device: " + str(master_data.iloc[row]["measurement_device"]))
		print("found in row: " + str(row))
		exit()
	else:
		sample_ID = sample_ID + "_" + str(master_data.iloc[row]["measurement_device"])
	if master_data.iloc[row]["gas"] not in gases:
		print("invalid gas: " + str(master_data.iloc[row]["gas"]))
		print("found in row: " + str(row))
		exit()
	else:
		sample_ID = sample_ID + "_" + str(master_data.iloc[row]["gas"])
		row_gases.update({row:str(master_data.iloc[row]["gas"])})
	row_ID.update({row:sample_ID})

print("will run the following sample IDs:")
print(row_ID)
#convert the sample ID column to be of string data type
master_data['Sample ID'] = master_data['Sample ID'].astype('object')
master_data["Use Data? (See Notes)"] = master_data["Use Data? (See Notes)"].astype('object')
print(master_data.dtypes)
for row, sample_ID in row_ID.items():
	master_data = analyze_slope(master_data, lgr_data,row,sample_ID,output_folder,r_2_values,row_gases[row],t_p_data)
	
#write the new master data file
master_data.to_excel(master_csv_path.replace('.xlsx','new.xlsx'), index = False)
exit()
