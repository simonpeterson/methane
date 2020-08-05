Input...script

>line 10: sample_id
#example: 2020_07-30_Vault-Lake_Bucket-T1 or BigTrail-Lake-Chamber-T1...{T2, T3, T4}...depending on lake and chamber used.
#remember to change "T1>T2>T3" each time. 

>line 15: sub_d (depth) {0,1,5,10} 
#unless, "collar", then comment out this following:
#collar_h = meancollar-3.53
#collar_h = np.mean([0])
*the first measuremnets you make 'FBX-CRRL-PT' will be bucket on 5/27/20
*but the next data uses *chamber* 'BTL-GSV-Collar' on 5/30/20

>line 19: (only change if you change "chamber type") {e.g. Chamber vs Bucket}
1. Bucket measurements use script section: 
# This section is optimized for the "bucket chambers" in saturated sediments. 
# snow density calculations are not required for summer (dont confuse winter with summer section)
2. Chamber measurements use script section:
# This section is optimized for the "large square chambers" rested on aluminum collar locations
*make sure to comment out the chamber, if using bucket, and vise versa

>line 25: #Atmospheric Pressure (inHg)
#https://www.wunderground.com/history/daily/us/ak/fairbanks/PAFA
#time of pressure reading: use 3pm for entire day.
P_inHg = 29.29 (uncomment) 

*just change this for each day.

FOR >line 38: 'date_time'...e.g. '6/11/2019 0:00' > (6/11/2019 2:30)
#comment out, we will use line 25 to derive our daily pressure/temp. 

>line 47: file, insert: gga_XXX file.txt > f0001 > f0002...etc. We only process "f000X files'

>line 55/56: start time, end time
#from excell doc

>line 74: CH4 or CO2
#I usually do CH4 first, then C02, remember to switch back to CH4.

Dont worry if you make mistake, we can always fix it. Good luck! 