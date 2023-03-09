# SENB-K_G_method
SENB - K and G (brittle polymer) method data analysis code for MATLAB

Code to determine K_IC and G_IC
In accordance to ISO 527-1 and ASTM D5045-14
Version 1.0
Written by Sammy He and Joe Terry
Department of Mechanical Engineering, Imperial College London
sch12@ic.ac.uk, jst114@ic.ac.uk
Copyright © 2023, Imperial College, London, All rights reserved
Last Updated: 09-03-2023


If you use this software for published work, use the following citation:

Terry, J.S., & He, S. (2022). SENB K and G (brittle polymer) method data analysis code for MATLAB.
Imperial College London, London, United Kingdom. [Software].  https://doi.org/


This code will find the values for K_IC, and G_IC (via energy and LEFM methods) for tests using multiple brittle plastic specimens, similar to that in ISO 527-1 and ASTM D5045 - 14

The code allows for a linear correction to the start of the load vs. extension data for compliance such that it passes through 0,0 and removes any initial load ramping errors.

The user can select the region to calculate the linear gradient for specimen compliance and to then calculate P5.

The code determines, Pmax, P5 and determines which to use. It then calculates energy under the curve for sample and compliance, and determines K_IC and G_IC (via energy and LEFM methods) for each specimen. 

It will apply the validity criteria, specify valid/invalid specimens and it will present the mean average values and standard deviation of K_IC and G_IC (via energy and LEFM methods) for valid specimens.

It is difficult to control initial crack length when razor tapping specimens, therefore the criteria specifying the length and variation of this length, can be slightly loosened if desired by the user.

INSTRUCTIONS

It is reccomended that you make a new copy of this matlab script for each type of specimen, so that it easy to re-run, should you so wish.

Specimen data files should be comma seperated value files (.csv). Load should be in N and displacement in mm, otherwise there will be discprepancies with the units in subsequent calculations. Specimen dimensions should be an excel spreadsheet (.xlsx) Dimension data should be inserted into a document following the format of 'Specimen Dimensions K G Template'. Dimensions should be in mm.

If you wish to add extra specimens to the dimensions excel sheet, you can copy and paste but do check the equation formulas are correct, especially those in row F which were attached to fixed cells e.g. $C$5, $C$9 etc.

To reduce probability of code errors, the compliance data should be put into the same folder as the specimen test data. File name for specimens should follow a system where the number at the end updates as the code iterates through the sample. e.g. 'Specimen_RawData_'. 

This will become Specimen_RawData_1, Specimen_RawData_2 etc. 
If preferred, the compliance (indentation) sample can follow this naming scheme, as the last specimen.

Often things can go wrong with testing (e.g. samples break before testing). For this code, the numbers of files must be consecutive. It is suggested to keep track of which samples correspond to the new naming system in a specimen_name file, if required for future reference.

In this matlab script, in the USER INPUT section:
- Write in a name for specimen type and add it's values for yield strength and Young's modulus from uniaxial tensile testing, these are required for validity checks and G_IC_LEFM.
- Insert file locations into the variables; directory_loc, indentdata_loc and dimension data.
- Input the start row of your data and the columns of load and extension from your specimen data files
- Input the threshold of load drop compared to maximum at which you would like to truncate test data, best to set this low e.g. 0.5 and increase if required for your particular dataset.
- Input the maxiumum load of compliance data to use during analysis
- Make any adjustment to the validity criteria of the initial crack length (e.g. extend from a/w = 0.45-0.55 or increase allowed variation from average initial crack length from 10%.

Then click run!

A figure will appear with the indentation data. Select the linear region that is representative of the sample compliance.

Afterwards, a figure will appear for each specimen, select two x-coordinates to define the range to determine the average linear gradient, this will be used to find the Load 5 line.
It is suggested to use this sensibly, select the linear section which is representative of the sample before the crack propagates. Some examples are shown in the help file.

Once the code has run, it will produce a written summary of the result for each specimen and mean values for valid specimens. It will also produce comparison figures for all specimens, displaying their validity.

If the user wishes to make the figure suitably sized for a report (8 cm wide), they may open the figure and then paste the three below lines into the command window.

set(gcf,'Units', 'centimeters','Position',[10 10 8 6])
set(gcf, [10 10 8 6]);
set(gca, 'XMinorTick','on','YMinorTick','on','Layer', 'top','Units', 'centimeters','Position',[0.9 0.9 6.8 4.8])
