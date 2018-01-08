
This is a brief read me for the latest version of 21CMMC (October 9th 2017).

This latest version contains the following:


1) Support for the simulated 21cm signal in both co-eval and light-cone formats
2) Inhomogeneous sub-grid recombinations
3) Ability to change the initial conditions on the fly within the MCMC
3) Ability to simultaneously constrain astrophysical and cosmological parameters
4) Line-of-sight redshift space distortions
5) Sampling of the 21cm global signal


21CMMC is a fast and efficient Python sampler of the semi-numerical reionisation 
simulation code 21cmFAST. It can recover constraints on astrophysical parameters 
from current or future 21 cm EoR experiments, accommodating a variety of EoR models, 
as well as priors on individual model parameters and the reionisation history. By 
studying the resulting impact on the EoR astrophysical constraints, 21CMMC can be 
used to optimise: (i) interferometer designs; (ii) foreground cleaning algorithms; 
(iii) observing strategies; (iv) alternate statistics characterising the 
21cm signal; and (v) synergies with other observational programs. 


For more details please see: http://adsabs.harvard.edu/abs/2015MNRAS.449.4246G
			     http://adsabs.harvard.edu/abs/2017MNRAS.472.2651G


Any questions, comments or feedback (including technical support and setting up) 
should be directed to me (Brad Greig) at brad.s.greig@gmail.com. Any issues will
aid me updating the code and user documentation!

If you use this code in any of your work, please cite the above papers. 


_______________________________________________________________________________


For any of the below topics, please refer to the corresponding documentation.


1) Setting up 21CMMC (“SettingUp21CMMC.txt”)

2) Preparing mock observations (“CreatingMockObservations.txt”)

3) Quick example setups of 21CMMC (“QuickStart.txt”) 

4) Brief descriptions of the available scripts (“ScriptDescriptions.txt”)


_______________________________________________________________________________


What is included in 21CMMC.

1) Example mock light-cone and co-eval 21cm PS and the 21cm global signal 
   in the directory “MockObs”

   These are to highlight naming conventions and file-formats etc.

2) Example light-cone, co-eval and 21cm global signal noise curves 
   in the directory “NoiseData”

   These are to highlight naming conventions and file-formats etc.

3) The neutral fraction PDFs for the Greig et al. (2017) analysis of the 
   z = 7.1 QSO (applicable as a prior).

