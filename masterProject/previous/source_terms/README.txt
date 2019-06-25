This code is a short calculation to get adequate source terms for CFD silmulation of a certain fluidised bed system. 

I have tried to be thorough in case someone without prior experience with python wants to use the code.

== Important to note before attempting to run the code:

This code is written in python 3, it will not work correctly if run using python 2. In addition, like most code it has dependencies that has to be installed first as well.

An easy way to install a python 3 distribution is to use a prepackaged distribution like anaconda or miniconda: https://conda.io/docs/user-guide/install/download.html

Dependencies:

numpy - for matlab-like math operations. Available in conda or pip repositories,quickfix: "conda install numpy" if using conda. otherwise google and follow the instuctions.

scipy - big library with a lot of functions for scientific computing. Needed for thermopy below. Can be installed similarily to numpy

thermopy - python3 library with thermodynamic tables for a wide variety of compounds. Has to be installed manually. easy install script but make sure you are inyour chosen python3 environment before running install. https://github.com/guillemborrell/thermopy

=======================================================================

How the scripts are meant to be used:

The main calculation is done by running the script "readcalc" it will read a casfile with inputs like gas yield and composition, take additional inputs specified in the top of the file and then print a csv with the source terms. It will also print a short summary to the console. 

If you want to provide a new case_file corresponding to another case than the chalmers gasifier you can generate that using the script "fuels-gas". to then use the new case_file you have to edit the corresponding read operation in readcalc.

Best of luck.

 
