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
![excel inputs](https://github.com/simonpeterson/methane/blob/master/readme_images/excel_inputs.PNG?raw=true)
https://github.com/simonpeterson/methane/blob/master/readme_images/excel_inputs.PNG

## 2. move the lgr data files into the input/lgr folder
The files **do not** need to be "cleaned", i.e. the pgp nonsense message that sometimes shows up at the end does not need to be deleted; the files are cleaned automatically.

## 3. temperature and pressure date
This will be updated later. Currently the pressure is pulled from a June 2019 file. This is something that can be changed at a later time; currently the program pulls the data from a
predetermined date, **not** the date of the measurement.

## 4. running the program
After all the desired samples have been entered, the program can be run. **Anaconda prompt** must be used to run the program; download this if you don't already have it. The program to be run is:

	> BTL_Collar_Monitring.py
	
there is also a .ipynb file of the same name; this is the Jupyter notebook file on which the program is based.

CD into the directory and run the program on anaconda prompt as shown:

![conda program run](https://github.com/simonpeterson/methane/blob/master/readme_images/conda_program_run.PNG?raw=true)