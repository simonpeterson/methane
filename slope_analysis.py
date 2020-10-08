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
import re

def analyze_slope(master_data, lgr_data,row,sample_ID,output_folder,r_values,gas_to_read,t_p_data):
	'''program that analyzes the slope and performs fitting'''
	#get the start and stop times from the file
	start_time = master_data.at[row, "start_time_(hh:mm:ss)"]
	stop_time = master_data.at[row, "stop_time_(hh:mm:ss)"]
	#the pressure data- if the start time is not available, pull from 2019 data
	#TOCHANGE- pull data from the actual date- is this a datetime object??
	P_Pa = t_p_data.loc[t_p_data['date_time'] == datetime.datetime.fromisoformat('2019-06-11'),'air_p_mean_Pa'].values
	
	#the measurement device used
	measurement_device = master_data.iloc[row]["measurement_device"]
	
	
	# #submersion depth in water or snow (cm)
	# #If multiple measurements are used in uneven surfaces, enter them in the mean field separated by comma
	# #when not submerged, (i.e. with collars) enter "0"
	sub_d = float(master_data.iloc[row]["submerged_depth(cm)"])
	# #Exposed height of chamber/bucket above surface (cm)
	# #If multiple measurements are used in uneven surfaces, enter them in the mean field separated by comma
	xh = float(master_data.iloc[row]["exposed_height(cm)"])-sub_d
	collar_h = float(master_data.iloc[row]["collar_height(cm)"])

	#TODO- add the pressure data to the file
	
	#lgr_datap = pd.read_csv(r'C:\Users\Simon\Documents\methane\June2019_T_P.csv', delimiter=',', parse_dates=[['date','time']])
	#P_Pa = lgr_datap.loc[lgr_datap['date_time'] == '6/11/2019 0:00','air_p_mean_Pa'].values
	
	#lgr_data.index = pd.DatetimeIndex(lgr_data.index)
	#lgr_data.index
	#lgr_data.keys()
	#lgr_data.head()
	
	#get the date to get the actual pressure dataS
	#master_data.iloc[row]["date_(yyyy-mm-dd)"]
	#the gas to read needs to be in the format:
	#[gas]d_ppm
	#while being right adjusted 20 spaces
	ts = ""
	ts          = lgr_data.between_time(start_time = "{}".format(start_time), end_time = "{}".format(stop_time))[('[' + gas_to_read + ']d_ppm').rjust(20)]
	temperature = lgr_data.between_time(start_time = "{}".format(start_time), end_time = "{}".format(stop_time))['              GasT_C']#.plot()
	#Plot the raw LGR data for the measurement window
	#change to check if ts is empty
	try:
		ts = xr.DataArray(ts, coords = [ts.index], dims = ['time'])
		ts.plot()
	except ImportError:
		print("the gives start time for data: " + sample_ID + "is not valid. will skip for now")
		master_data.at[row, "program_run?"] = "n"
		return master_data
	new_ts = ts.dropna('time')
	if new_ts.time.size == 0:
		print("the gives start time for data: " + sample_ID + "is not valid. will skip for now")
		return master_data
	temperature_mean = xr.DataArray(temperature, coords = [temperature.index], dims = ['time']).mean().data
	temperature_error = temperature.std() + 1
	#lgr_data.head()
	type(temperature_mean)

	# In[75]:


	# #This section is optimized for the "bucket chambers" in the snow. 
	# # #A snow density is assumed from literature values, and the snow water equivalent is calculated to 
	# # #correct for the snow's porosity in the chamber
	# # #Conditions for snow or sediment height estimation
	if xh <= 5:
	    xh_error = xh*0.5
	else :
	     xh_error = xh*0.2
	# #Volume of Chamber (as a truncated cone) (cm^3)
	# vol = pi*h*(r_bottom**2 + r_top**2 + r_bottom*r_top)/3

	#Normal Height (height from center of bottom base to apex of the imaginary cone) (cm)
	# n_height = h + ((h*r_top)/(r_bottom - r_top))

	#Volume of Normal cone (if chamber were not truncated, i.e. bucket-shaped)(cm^3)
	# cone_volume = (1/3)*pi*(r_bottom**2)*n_height

	#Volume of the imaginary cone (the volume that is removed/truncated to form the bucket shape)(cm^3)
	# imagine_volume = (1/3)*pi*(r_top**2)*(n_height - h)

	#Hypotenuse of Normal cone (cm)
	# hypo_cone = sqrt((r_bottom**2)+(n_height**2))

	#angle of normal cone apothems (radians)
	# angle = (asin(r_bottom/hypo_cone))*2

	#Normal height of Snow Cone (cm)
	# n_height_snow = (n_height - h) + (h - xh)
	# n_height_snow_error = xh_error

	#exterior angle *between snow level and the side of bucket (downward) (degrees)
	# ex_angle = 180 - ((angle*(180/pi))/2) - 90

	# # #Hypotenuse of Snow cone (cm)
	# hypo_snow = (n_height_snow/cos(angle/2))
	# hypo_snow_error = n_height_snow_error*abs((cos(angle/2)))

	# # #Radius of top of snow (cm)
	# r_snow = sqrt((hypo_snow**2)-(n_height_snow**2))
	# r_snow_error = 2*(hypo_snow_error/abs(hypo_snow))*abs(r_snow) + 2*(n_height_snow_error/abs(n_height_snow))*abs(r_snow)

	# # #Volume of snow cone (cm^3)
	# snow_cone_volume = (1/3)*pi*(r_snow**2)*(n_height_snow)
	# snow_cone_volume_error = (sqrt((r_snow_error/r_snow)**2 + (n_height_snow_error/n_height_snow)**2))*abs(pi/3)

	# # #Uncorrected Volume of Snow in Chamber (non-porous snow) (cm^3)
	# un_snow_vol = snow_cone_volume - imagine_volume
	# un_snow_vol_error = snow_cone_volume_error

	# # #Snow bulk density (g/cm^3) Sturm et al. 2010 Journal of Hydrometeorology
	# # #Mean = 0.217 g/cm^3 Taken from 1541 observations in Alaskan/Canadian Taiga
	# # #Std = 0.056 g/cm^3 
	# pb = 0.217
	# pb_error = 0.056

	# # #Snow water equivalent (cm)
	# swe = (h - xh)*(pb/1)
	# swe_error = sqrt((xh_error/xh)**2 + (pb_error/pb)**2)

	# # #Corrected volume of snow (cylinder since "top" base of bucket is pushed into snow)(cm^3)
	# corrected_snow_vol = pi*(r_top**2)*swe
	# corrected_snow_vol_error = swe_error *abs(pi*(r_top**2))

	# # #Total volume(cm^3)
	# total_vol = vol -un_snow_vol + lid_vol + t_vol + LGR_vol
	# total_vol_error = un_snow_vol_error + lid_vol_error + t_vol_error + LGR_vol_error

	# # #Total volume corrected (accounting for swe-based snow volume)(cm^3)
	# total_vol_corr = vol - corrected_snow_vol + lid_vol + t_vol + LGR_vol
	# total_vol_corr_error = corrected_snow_vol_error + lid_vol_error+ t_vol_error + LGR_vol_error

	# # #Convert total corrected volume to liters (L)
	# V = total_vol_corr / 1000
	# V_error = total_vol_corr_error / 1000

	# # #Surface area of chamber opening (m^2)
	# area = pi*(r_top**2)/10000


	# In[46]:

	
	# This section is optimized for the "bucket chambers" in saturated sediments. 
	# snow density calculations are not required 
	print("calculating size")
	
	
	#if the container is a bucket, define variables that are used for both sediment and snow measurements
	if "bucket" in master_data.iloc[row]["measurement_device"]:
		# #lid volume (cm^3)
		lid_vol = 975
		lid_vol_error = lid_vol* 0.01
		#dimensions of the bucket
		h = 34.5
		diam_top = 25.8
		r_top = diam_top/2
		# #tubing dimensions (cm, cm^3 for t_vol)
		t_length = 914.4
		t_length_error = t_length * 0.05
		t_id = 0.3175
		t_vol = pi*((t_id/2)**2)*t_length
		t_vol_error = t_length_error * abs(t_id)
		# #LGR cell volume (cm^3)
		LGR_vol = 335
		LGR_vol_error = LGR_vol*0.01
		# #"bottom" here is actually the top of the bucket. Because when modelling a truncated cone, it's 
		# #easier to think of the bottom as the base with the larger radius
		diam_bottom = 29
		r_bottom = diam_bottom/2
		if master_data.iloc[row]["measurement_device"] == "bucket_sediment":
			print("using bucket in sediment calculation")

			#Volume of Chamber (as a truncated cone) (cm^3)
			vol = pi*h*(r_bottom**2 + r_top**2 + r_bottom*r_top)/3

			#Normal Height (height from center of bottom base to apex of the imaginary cone) (cm)
			n_height = h + ((h*r_top)/(r_bottom - r_top))

			#Volume of Normal cone (if chamber were not truncated, i.e. bucket-shaped)(cm^3)
			cone_volume = (1/3)*pi*(r_bottom**2)*n_height

			#Volume of the imaginary cone (the volume that is removed/truncated to form the bucket shape)(cm^3)
			imagine_volume = (1/3)*pi*(r_top**2)*(n_height - h)

			#Hypotenuse of Normal cone (cm)
			hypo_cone = sqrt((r_bottom**2)+(n_height**2))

			#angle of normal cone apothems (radians)
			angle = (asin(r_bottom/hypo_cone))*2

			#Normal height of sediment Cone (cm)
			n_height_sed = (n_height - h) + (h - xh)
			n_height_sed_error = xh_error

			#exterior angle *between sediment level and the side of bucket (downward) (degrees)
			ex_angle = 180 - ((angle*(180/pi))/2) - 90

			#Hypotenuse of sediment cone (cm)
			hypo_sed = (n_height_sed/cos(angle/2))
			hypo_sed_error = n_height_sed_error*abs((cos(angle/2)))

			#Radius of top of sediment (cm)
			r_sed = sqrt((hypo_sed**2)-(n_height_sed**2))
			r_sed_error = 2*(hypo_sed_error/abs(hypo_sed))*abs(r_sed) + 2*(n_height_sed_error/abs(n_height_sed))*abs(r_sed)

			#Volume of sediment cone (cm^3)
			sed_cone_volume = (1/3)*pi*(r_sed**2)*(n_height_sed)
			sed_cone_volume_error = (sqrt((r_sed_error/r_sed)**2 + (n_height_sed_error/n_height_sed)**2))*abs(pi/3)

			#Uncorrected Volume of sediment in Chamber (non-porous snow) (cm^3)
			un_sed_vol = sed_cone_volume - imagine_volume
			un_sed_vol_error = sed_cone_volume_error

			#Snow bulk density (g/cm^3) Sturm et al. 2010 Journal of Hydrometeorology
			#Mean = 0.217 g/cm^3 Taken from 1541 observations in Alaskan/Canadian Taiga
			#Std = 0.056 g/cm^3
			pb = 0.217
			pb_error = 0.056

			#Snow water equivalent (cm)
			swe = (h - xh)*(pb/1)
			swe_error = sqrt((xh_error/xh)**2 + (pb_error/pb)**2)

			#Corrected volume of snow (cylinder since "top" base of bucket is pushed into snow)(cm^3)
			corrected_snow_vol = pi*(r_top**2)*swe
			corrected_snow_vol_error = swe_error *abs(pi*(r_top**2))

			#Total volume(cm^3)
			total_vol = vol -un_sed_vol + lid_vol + t_vol + LGR_vol
			total_vol_error = un_sed_vol_error + lid_vol_error + t_vol_error + LGR_vol_error

			#Total volume corrected (accounting for swe-based snow volume)(cm^3)
			total_vol_corr = vol - corrected_snow_vol + lid_vol + t_vol + LGR_vol
			total_vol_corr_error = corrected_snow_vol_error + lid_vol_error+ t_vol_error + LGR_vol_error

			#Convert total corrected volume to liters (L)
			V = total_vol / 1000
			V_error = total_vol_error / 1000

			#Surface area of chamber opening (m^2)
			area = pi*(r_top**2)/10000
		elif master_data.iloc[row]["measurement_device"] == "bucket_snow":
			print("running snow bucket calculation")
			#Volume of Chamber (as a truncated cone) (cm^3)
			vol = pi*h*(r_bottom**2 + r_top**2 + r_bottom*r_top)/3

			#Normal Height (height from center of bottom base to apex of the imaginary cone) (cm)
			n_height = h + ((h*r_top)/(r_bottom - r_top))

			#Volume of Normal cone (if chamber were not truncated, i.e. bucket-shaped)(cm^3)
			cone_volume = (1/3)*pi*(r_bottom**2)*n_height

			#Volume of the imaginary cone (the volume that is removed/truncated to form the bucket shape)(cm^3)
			imagine_volume = (1/3)*pi*(r_top**2)*(n_height - h)

			#Hypotenuse of Normal cone (cm)
			hypo_cone = sqrt((r_bottom**2)+(n_height**2))

			#angle of normal cone apothems (radians)
			angle = (asin(r_bottom/hypo_cone))*2

			#Normal height of Snow Cone (cm)
			n_height_snow = (n_height - h) + (h - xh)
			n_height_snow_error = xh_error

			#exterior angle *between snow level and the side of bucket (downward) (degrees)
			ex_angle = 180 - ((angle*(180/pi))/2) - 90

			#Hypotenuse of Snow cone (cm)
			hypo_snow = (n_height_snow/cos(angle/2))
			hypo_snow_error = n_height_snow_error*abs((cos(angle/2)))

			#Radius of top of snow (cm)
			r_snow = sqrt((hypo_snow**2)-(n_height_snow**2))
			r_snow_error = 2*(hypo_snow_error/abs(hypo_snow))*abs(r_snow) + 2*(n_height_snow_error/abs(n_height_snow))*abs(r_snow)

			#Volume of snow cone (cm^3)
			snow_cone_volume = (1/3)*pi*(r_snow**2)*(n_height_snow)
			snow_cone_volume_error = (sqrt((r_snow_error/r_snow)**2 + (n_height_snow_error/n_height_snow)**2))*abs(pi/3)

			#Uncorrected Volume of Snow in Chamber (non-porous snow) (cm^3)
			un_snow_vol = snow_cone_volume - imagine_volume
			un_snow_vol_error = snow_cone_volume_error

			#Snow bulk density (g/cm^3) Sturm et al. 2010 Journal of Hydrometeorology
			#Mean = 0.217 g/cm^3 Taken from 1541 observations in Alaskan/Canadian Taiga
			#Std = 0.056 g/cm^3 
			pb = 0.217
			pb_error = 0.056

			#Snow water equivalent (cm)
			swe = (h - xh)*(pb/1)
			swe_error = sqrt((xh_error/xh)**2 + (pb_error/pb)**2)

			#Corrected volume of snow (cylinder since "top" base of bucket is pushed into snow)(cm^3)
			corrected_snow_vol = pi*(r_top**2)*swe
			corrected_snow_vol_error = swe_error *abs(pi*(r_top**2))

			#Total volume(cm^3)
			total_vol = vol -un_snow_vol + lid_vol + t_vol + LGR_vol
			total_vol_error = un_snow_vol_error + lid_vol_error + t_vol_error + LGR_vol_error

			#Total volume corrected (accounting for swe-based snow volume)(cm^3)
			total_vol_corr = vol - corrected_snow_vol + lid_vol + t_vol + LGR_vol
			total_vol_corr_error = corrected_snow_vol_error + lid_vol_error+ t_vol_error + LGR_vol_error

			#Convert total corrected volume to liters (L)
			V = total_vol_corr / 1000
			V_error = total_vol_corr_error / 1000

			#Surface area of chamber opening (m^2)
			area = pi*(r_top**2)/10000
		else:
			print("invalid measurement container" + master_data.iloc[row]["measurement_device"] + " in row: " + str(row))
			master_data.at[row, "program_run?"] = "n"
			return master_data
	elif master_data.iloc[row]["measurement_device"] == "chamber":
		# This section is optimized for the "large square chambers" rested on aluminum collar locations. 
		# density calculations are not required 
		print("running chamber")
		#Chamber dimensions for large clear chamber (units of cm)
		h = 106.5 - float(sub_d)
		w = 65.7
		d = 65.7
		vol = h*w*d
		vol_error = vol * 0.05
		xh = h + collar_h

		#Collar dimensions (units of cm)
		collar_w = 70.7
		collar_d = 70.7
		collar_vol = collar_h*collar_w*collar_d
		collar_vol_error = collar_vol * 0.05

		#tubing dimensions (cm, cm^3 for t_vol)
		t_length = 914.4
		t_length_error = t_length * 0.05
		t_id = 0.3175
		t_vol = pi*((t_id/2)**2)*t_length
		t_vol_error = t_length_error * abs(t_id)

		#LGR cell volume (cm^3)
		LGR_vol = 335
		LGR_vol_error = LGR_vol*0.01

		#Snow bulk density (g/cm^3) Sturm et al. 2010 Journal of Hydrometeorology
		#Mean = 0.217 g/cm^3 Taken from 1541 observations in Alaskan/Canadian Taiga
		#Std = 0.056 g/cm^3
		pb = 0.217
		pb_error = 0.056

		#Snow water equivalent (cm)
		swe = (h - xh)*(pb/1)
		swe_error = sqrt((xh_error/xh)**2 + (pb_error/pb)**2)
		
		
		#ACHTUNG!!!!
		#this was commented out by Simon Peterson on 27.09.20, as it does not make sense for snow errror 
		#calculations to be in the chamber section
		#Corrected volume of snow (cylinder since "top" base of bucket is pushed into snow)(cm^3)
		#corrected_snow_vol = pi*(r_top**2)*swe
		#corrected_snow_vol_error = swe_error *abs(pi*(r_top**2))

		#Total volume(cm^3)
		total_vol = vol + t_vol + LGR_vol + collar_vol
		total_vol_error = vol_error + t_vol_error + LGR_vol_error + collar_vol_error

		#Convert total corrected volume to liters (L)
		V = total_vol / 1000
		V_error = total_vol_error / 1000

		#Surface area of collar opening (m^2)
		area = (collar_w*collar_d)/10000
	else:
		print("invalid measurement container" + master_data.iloc[row]["measurement_device"] + " in row: " + str(row))
		master_data.at[row, "program_run?"] = "n"
		return master_data
	# As a sanity check, print the volume of the chamber/bucket and it's error in liters
	print('volume = ', V, ' ±',  V_error,  ' liters')
	
	# Goal: 
	# We find sections in the observation window which are
	# 1. at least 50 seconds long and no more than 210 (50 < n > 210)
	# 2. the linear slope of the section should fit the data with R2 0.985. WHich means 98% of the variation in data can be explained by our model. 
	# 
	# Algorithm:
	# We start with a section of first 210 points from the section window.
	# We fit a line to the section and compute the errors as the squareroot of the diagnols of the covariance matrix between the model and data.
	# We also compute R2 between the slope and the data and if R2 is greater than 0.985, we stop, since we got the largest continuous section which fulfils the requirement.
	# i.e. our best case scenario! If we don't get a valid slope value, we shift the section by 1. i.e. from the second element to 211th element in the section window and repeat the process.
	# If we don't get a valid slope even after shifting through the entire length of the section window, we reduce the section length to 209 and start sliding. 
	# So on until we reduce the section length to its lowest threshold, i.e. 50. 
	# 
	# If we still don't get a valid value for the smallest section, we start smoothing the sections by 5, 10, and ultimately 15-point moving windows. 
	#CHANGES BY SIMON PETERSON 09.09.20
	# Make the R2 value 
	# In[]
	# Set data quality thresholds
	min_section_length_original = 45
	max_section_length = 210
	# r_2_threshold = 0.985
	#starting smoothing window
	smoothing_window = 0
	########################################################################
	section_length   = np.arange(min_section_length_original,max_section_length+1)[::-1]
	for r_2_value in r_values:
		max_r_2 = []
		print("running with r^2 of: " + str(r_2_value))
		valid_section_length, slope, a, R_squared, max_r_squared = brain(section_length, ts,r_2_value )
		max_r_2.append(max_r_squared)
		if slope:
			print('valid section length = %.3f' %valid_section_length)
			print('Smoothing_window = %.3f' %smoothing_window)
	#                 print('Slope = %.4f' %slope)
	#                 print('Temperature = %.3f' %temperature_mean)
	#                 print('Section start timestamp = ' +str(a[0].time.data))
			break
		else:
			print("running with smoothing window of 5")
			smoothing_window  = 5
			min_section_length = min_section_length_original*2
			section_length   = np.arange(min_section_length,max_section_length+1)[::-1]
			ts_smooth = ts.rolling(time = smoothing_window, center = True).mean().dropna(dim='time')
			valid_section_length, slope, a, R_squared, max_r_sqaured = brain(section_length, ts,r_2_value)
			max_r_2.append(max_r_squared)
			if slope:
				print('valid section length = %.3f' %valid_section_length)
				print('Smoothing_window = %.3f' %smoothing_window)
		#                 print('Slope = %.4f' %slope)
		#                 print('Temperature = %.3f' %temperature_mean)
		#                 print('Section start timestamp = ' +str(a[0].time.data))    
				break
			else:
				print("running with smoothing window of 10")
				smoothing_window  = 10
				min_section_length = int(min_section_length_original*2)
				section_length   = np.arange(min_section_length,max_section_length+1)[::-1]
				ts_smooth = ts.rolling(time = smoothing_window, center = True).mean().dropna(dim='time')
				valid_section_length, slope, a, R_squared, max_r_squared = brain(section_length, ts,r_2_value )
				max_r_2.append(max_r_squared)
				if slope:
					print('valid section length = %.3f' %valid_section_length)
					print('Smoothing_window = %.3f' %smoothing_window)
		#                 print('Slope = %.4f' %slope)
		#                 print('Temperature = %.3f' %temperature_mean)
		#                 print('Section start timestamp = ' +str(a[0].time.data))        
					break
				else:
					smoothing_window  = 15
					min_section_length = min_section_length_original*3
					section_length   = np.arange(min_section_length,max_section_length+1)[::-1]
					ts_smooth = ts.rolling(time = smoothing_window, center = True).mean().dropna(dim='time')
					valid_section_length, slope, a, R_squared, max_r_squared = brain(section_length, ts,r_2_value )
					max_r_2.append(max_r_squared)
					if slope:
						print('valid section length = %.3f' %valid_section_length)
						print('Smoothing_window = %.3f' %smoothing_window)
		#                 print('Slope = %.4f' %slope)
		#                 print('Temperature = %.3f' %temperature_mean)
		#                 print('Section start timestamp = ' +str(a[0].time.data))
						break
					elif r_2_value == min(r_values):
						print('Didn\'t work with lowest R_2 threshold value! Baaad data!!! will still')
						master_data.at[row, "program_run?"] = "y"
						master_data.at[row, "Use Data? (See Notes)"] = "rejected"
						master_data.at[row, "R_value_used"] = max(max_r_2)
						return master_data
	print('valid section length = %d' %valid_section_length)
	print('Smoothing_window = %d' %smoothing_window)
	print('Slope = %.5f ppm/s' %slope)
	print('R^2 = %.4f ' %R_squared[0].data)
	print('R^2 = %.4f ' %R_squared[0].data)
	print('Temperature = %.3f deg. C' %temperature_mean)
	print('Section start timestamp = ' +str(a[0].time.data))

	#Convert the section start time from datetime stamp to actual index of the original section window.
	index = np.where(ts.time ==a[0].time)[0][0]
	slope, slope_error, ts_section, y, y_hat = get_slope_error(ts[index: index+valid_section_length], plot = True)
	#print('Slope error = %.6f ppm/s' %slope_error)

	#Calculate Flux and its +/- uncertainty (micromoles m^-2 s^-1)
	#Gas Constant L atm mol^-1 K^-1
	#pressure in atmospheres (atm)
	P = 9.86923e-6 * P_Pa
	P_error = P * 0.01
	R = 0.082057338
	moles = (P*V)/(R*(temperature_mean+273.15))
	moles_error = sqrt((P_error/P)**2 + (V_error/V)**2 + (temperature_error/(temperature_mean+273.15)**2))
	flux = slope*moles/area
	flux_error = (sqrt(slope_error/abs(slope))**2 + (moles_error/moles)**2)
	print('Flux = %.3f ' %flux + '± %.3f ' %flux_error)
	#print('Flux Error = ± %.3f ' %flux_error)
	print('Units = micromol CH4 m^-2 s^-1')

	#make the Figure
	figure, ax = plt.subplots(figsize = (10,8))
	ax.plot(ts_section.time, y)
	ax.plot(ts_section.time, y_hat)
	plt.ylabel('PPM CH4_$', fontsize = 12)
	plt.xlabel('Time',fontsize = 12)
	legend = plt.legend(('Raw Data', 'Linear Fit'), title = sample_ID, loc='upper left', fontsize = 12, shadow=True)
	plt.setp(legend.get_title(),fontsize='large', fontweight = 'bold')
	plt.text(0.02,0.7, 'Flux = %.3f ± %.3f \nmicromol $CH_4$ $m^{-2}$ $s^{-1}$ \nR$^2$ = %.3f' %(flux, flux_error,R_squared[0].data), fontsize = 13, transform=ax.transAxes)

	# #save the figure. 
	# #########!!!INPUT!!!########!!!INPUT!!!####!!!INPUT!!!##################
	# #########!!!INPUT!!!########!!!INPUT!!!####!!!INPUT!!!##################
	print("saving output file as: \n" + str(output_folder + "\\" + sample_ID + '.png'))
	plt.savefig(str(output_folder + "\\" + sample_ID + '.png'),dpi = 300, bbox_inches = 'tight')
	# ########################################################################
	# print('xh = %.1f' %xh)
	# Out put the data from the computation. 
	# With each program iteration, the results from the computation are saved into the space-delimited
	# .txt file specified below (i.e. "July2019FluxData.txt")
	# output variables to save into the output log file
	# output_data = [sample_id,temperature_mean,P[0],area,V,xh,valid_section_length,smoothing_window,slope,slope_error,R_squared[0].data,str(a[0].time.data),flux[0],flux_error[0]]
	#########!!!INPUT!!!########!!!INPUT!!!####!!!INPUT!!!##################
	#########!!!INPUT!!!########!!!INPUT!!!####!!!INPUT!!!##################

	# ########################################################################  
	# print(output_data)
	# use this space delimited text (staring with "Sample_ID") to fill in the header line once time after the file is generated
	# for the first time.
	# NOTE: The header labels must be in the same order as the data in the "output_data" array:
	# Sample_ID Temp_C P_atm Area_m2 ChamberVol_L Chamber_height_cm valid_section_seconds smoothing_seconds slope_ppm_-s slope_error Rsquared Start_time flux_umol_m-2_s flux_error


	# In[67]:3


	# Out+ put the data from the computation. 
	# With each program iteration, the results from the computation are saved into the space-delimited
	# .txt file specified below (i.e. "July2019FluxData.txt")

	# output variables to save into the output log file
	output_data_headers = ["sample_ID", "Temperature","Pressure","area", "volume", "valid_section_length","smoothing_window","slope", "slope error",
	"R^2","time","flux","flux error"]
	output_data = [sample_ID.replace("꞉",":"),temperature_mean,P[0],area,V,xh,valid_section_length,smoothing_window,slope,slope_error,
				   R_squared[0].data,str(a[0].time.data),flux[0],flux_error[0]]	
	#updating the master spreadsheet
	#update the sample ID
	master_data.at[row, "Sample ID"] = sample_ID
	#the pressure used in the measurements, measured in pascals??
	master_data.at[row, "air_Pa"] = P[0]
	#the r squared value used
	master_data.at[row, "R_value_used"] = R_squared[0]
	#the flux found- need to know which gas we are measuring
	if gas_to_read == "CH4":
		master_data.at[row, "CH4 flux μmol m^-2 s^-1"] = flux[0]
		master_data.at[row, "CH4 flux ± uncertainty"] = flux_error[0]
	elif gas_to_read == "CO2":
		master_data.at[row, "CO2 Flux μmol m^-2 s^-1"] = flux[0]
		master_data.at[row, "CO2 flux ± uncertainty"] = flux_error[0]
	else:
		print(gas_to_read + " not able to be run")
		return master_data
	master_data.at[row, "program_run?"] = "y"
	#
	#########!!!INPUT!!!########!!!INPUT!!!####!!!INPUT!!!##################
	#########!!!INPUT!!!########!!!INPUT!!!####!!!INPUT!!!##################
	with open(os.path.join(output_folder,(sample_ID + ".txt")),'w') as f: 	
		np.savetxt(f,output_data_headers,fmt = "%s,",newline=' ')
		f.write('\n')
		np.savetxt(f, output_data, fmt = "%s,",delimiter=',', newline=' ')
	########################################################################
		
	#edit the master data sheet with the data from the 
	return master_data
def compute_r2(ts_section, plot = False):
     x     = np.arange(len(ts_section))
     #removes time metadata for simplicity
     y     = np.array(ts_section)
     model = np.polyfit(x,y,1)
     slope = model[0]
     intercept = model[1]
     y_hat = slope*x + intercept # OR np.polyval(model, x)
     correlation_coefficient = np.corrcoef(y,y_hat)[0,1]
     r_square = (correlation_coefficient)**2
     if plot:
         plt.plot(y)
         plt.plot(y_hat)

     return r_square, slope
def brain(section_length, ts,r_2_threshold):
    max_r_2_values = []
    for section_length in section_length:
        # print(section_length)
        r_2 = []
        slope = []
        start = []
        # end   = []
        # section_length = 210
        for i in range(len(ts)-section_length):
            ts_section = ts[i:i + section_length]
            tmp_start  = ts_section.time[0]
            tmp_end    = ts_section.time[-1]
            tmp_r2, tmp_slope = compute_r2(ts_section, plot=False)
            r_2.append(tmp_r2)
            slope.append(tmp_slope)
            start.append(tmp_start)
            # end.append(tmp_end)
        r_2 = xr.DataArray(r_2, coords = [ts[:-int(section_length)].time], dims = ['time'])
        slope = xr.DataArray(slope, coords = [ts[:-int(section_length)].time], dims = ['time'])
        # end = xr.DataArray(end, coords = [ts[:-int(section_length)].time], dims = ['time'])
        # start = xr.DataArray(start, coords = [ts[:-int(section_length)].time], dims = ['time'])
        # just select those r2s that qualify the threshold, drop the rest to get time stamp

        # c = end.where(r_2>=r_2_threshold, drop=True)
        #
        #all of the data windows that fit the r_2 threshhold
        a = list(slope.where(r_2>=r_2_threshold, drop=True))
		#all the R_squared values that worked
        R_squared = list(r_2.where(r_2>=r_2_threshold, drop=True))
		#the maximum r_squared value
		#doesn't work if the array size is 0
        if (len(ts)-section_length) <= 0:
            max_r_sqaured = 0
        else:
            max_r_squared_index = r_2.argmax()
            max_r_squared = r_2[max_r_squared_index]
            max_r_2_values.append(max_r_squared)
        result = []
        valid_length = []
        if a:
            # print(a)
            result = a[0].data
            valid_length = section_length
            # print('starting time stamp' + a[0].time)
            break

        # ts.sel(time=slice(b[0],c[0]))
        # ts.plot()
        # ts.loc[a.time.values].plot()
        #
        # np.where(ts.index==a.time[0])
    return valid_length, result, a, R_squared, max(max_r_2_values)

def get_slope_error(ts_section, plot = True):
     x     = np.arange(len(ts_section))
     #removes time metadata for simplicity
     y     = np.array(ts_section)
     model, M = np.polyfit(x,y,1, cov = True)
     slope = model[0]
     slope_error = np.sqrt(M[0][0])
     intercept = model[1]
     y_hat = slope*x + intercept # OR np.polyval(model, x)

     if plot:
            print('Slope error = %.6f ppm/s' %slope_error)
#             plt.figure(figsize = (10,8))
#             plt.plot(ts_section.time, y)
#             plt.plot(ts_section.time, y_hat)
     return slope, slope_error, ts_section, y, y_hat
