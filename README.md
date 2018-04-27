# -RRC-Heng-Tsai-2016-ApJ
TEA compendium for Heng &amp; Tsai(2016), "Analytical Models of Exoplanetary Atmospheres III: Gaseous C–H–O–N Chemistry with Nine Molecules"

This is the TEA compendium for the paper, Heng & Tsai (2016) "Analytical Models of Exoplanetary Atmospheres. III. Gaseous C-H-O-N Chemistry with 9 Molecules". The paper is submitted to the ApJ in March 17, 2016.
The paper can be found via http://iopscience.iop.org/article/10.3847/0004-637X/829/2/104/pdf or doi:10.3847/0004-637X/829/2/104

The purpose of this code is to validate our analytical formula against the open-source code TEA ([Blecic et al. (2015)](https://arxiv.org/abs/1505.06392)), which numerically
calculates equilibrium chemistry by Gibbs free energy minimization. We demonstrate that the analytical formula is accurate at the ∼ 1% level for temperatures from 500-3000 K.

The TEA version used to produce the results from this paper is located at:
https://github.com/dzesmin/TEA

The files in this directory are named according to the figures from the paper.
There are seven .atm files as the input files for TEA, seven .tea pre-run output files from TEA, and three python files to make Figure 1. 2. 3. in the paper. The following instructions show how to run the code and reproduce the figures.


1. Each input and output files are named according to the figure numbers. For example, "f1a.atm" and "f1a.tea" are the input and output files of TEA for Figure 1a (the top panel). For a quick plot, skip to Step 3. to generate the plot from the pre-run output files. Step 3. is to re-produce the output files. 
 
2. Run "runatm.py" with the corresponding atm file, as the standard preocedure of performing multiple T,P calculation for TEA. Users are referred to the manual of TEA "TEA-UserManual.pdf". Copy the output file (.tea) back into the current directory.  

3. With the TEA output file in the same directory as the python script "heng_tsai_fX.py", where X is again the number of the figures, open the python script and change the first line 

plot_fig = '1a' 

according to the current figure to plot, e.g. plot_fig = '1c' for the bottom plot in Figure 1.
Finally, run the python script to generate the plot. 

More details about using TEA can be found on https://github.com/dzesmin/TEA