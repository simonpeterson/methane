# greenhouse gas flux analysis program
This edited version of code originally written by JPL allows for fast and efficient analysis of greenhouse gas flux data. Written by Simon Peterson for the Katey Anthony lab at the University of
Alaska Fairbanks, started on 08/24/2020
## Usage instructions:
There is a master data excel file called **simon_masters.xlsx** where the data is entered, located in the data subfolder. For the program to be run, the following data must be present: 

- start time
- stop time
- location
- date
- measurement device
- collar
- collar height
- submerged depth
- exposed height

These values must be entered for **each desired measurement time**
## 1. enter the measurement time(s) into the master excel file
For example, a measurement with the following data:

- **start time:** 			17:12:20
- **stop time:**  			17:15:30
- **location**    			vault lake
- **date**        			2020-07-30
- **measurement device** 	bucket
- **collar**             	yes
- **collar height**         0 cm
- **submerged depth**       0 cm
- **exposed height**        34.5 cm

would look like this:
![excel inputs](https://github.com/simonpeterson/methane/blob/master/readme_images/excel_inputs.PNG?raw=true)
The `<program_run?>`
 column does not need to have an "n"; anything that isn't "y" will cause the program to process the data
## 2. Move the lgr Data Files into the input/lgr Folder
The files **do not** need to be "cleaned", i.e. the pgp nonsense message that sometimes shows up at the end does not need to be deleted; the files are cleaned automatically.

## 3. Temperature and Pressure Data
This will be updated later. Currently the pressure is pulled from a June 2019 file. This is something that can be changed at a later time; currently the program pulls the data from a
predetermined date, **not** the date of the measurement.

## 4. Running The Program
After all the desired samples have been entered, the program can be run. **Anaconda prompt** must be used to run the program; download this if you don't already have it. The program to be run is:

	> BTL_Collar_Monitring.py
	
there is also a .ipynb file of the same name; this is the Jupyter notebook file on which the program is based.

CD into the directory and run the program on anaconda prompt as shown:

![conda program run](https://github.com/simonpeterson/methane/blob/master/readme_images/conda_program_run.PNG?raw=true)

The program will then run. There will be many messages flying across the screen; many of these are currently for debugging purposes, and can be ignored.

## 5. Outputs
There are three main "classes" of output files:

1. A new excel file is created which is simply the name of the old file + "_new". This file can be reviewed for accuracy then saved as the master data file. It is located in the **data** subfolder. Currently, this file is "simon_masters_new.xlsx"
2. A plot of the fit. This is done for each gas for each sample time given in the excel file. These are found in the **outputs** subfolder.
3. A text output of the data. This is also put into the "outputs" file. 

The image and the text output data files are named after the sample ID, which takes the following format:

	> yyyy-mm-dd_HHh-mmM-ssS_location_collection-device_gas

where:
- **yyyy-mm-dd** is the date of the sample
- **hhH-mmM-ssS** is the start time of the sample
- **location** is the location where the sample was collected
- **collection-device** is the device used to gather the data (bucket, chamber, etc.)
- **gas** is the gas that was measured.

an example of the output files is shown:
![example outputs](https://github.com/simonpeterson/methane/blob/master/readme_images/example_outputs.PNG?raw=true)

