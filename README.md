# Assigment-4
Project 4 in FYS4150

This project is using the two dimensional Ising model and the Metropolis algorithm to calculate thermodynamic variables. Most variables are generated in the c++ part, while some, like the variance and energy distribution, is calculated in MatLab.

The c++ code uses armadillo, so remember to change the path in the pro-file to where your armadillo library are located!
MPI has been used to parallellize the code, so you will also have to change some settings for this to work on your computer. In this project, eight processors were used.
When using the MatLab code, remember to change the path of where MatLab looks for file. This is done by "addpath('your path here at the top of each MatLab script. The files you need are already there, MatLab only needs to find them.
The plots themselfes are saved as .jpg files in the article folder to be used in LaTeX. You will find more plots in the folder, than in the pdf, because some plots never ended up being used.
