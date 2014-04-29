CSE-6730-Project-2
==================
The primary files in this package are as follows
MKS - module containing methods to execute MKS
AbaqusGen - module with functions to create Abaqus input files
Plot - module for plotting and analyzing parts of microstructure
MSf - function to create microstructure function from input microstructure

TestAbaqusGen - function used to set parameters for Abaqus input file to generate
predictSteadyStateMKs - script to run MKS on different calibrations and responses
	This script can be run and the parameters adjusted to calculate a desired response

All of these packages require numpy and scipy to be installed to run properly.
Abaqus output files and respective inputs may be found in the calibration and validation modeling folders


