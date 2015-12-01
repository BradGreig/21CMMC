This is a (likely, highly incomplete) read me file for the first (latest) release of 21CMMC!

21CMMC is a fast and efficient Python sampler of the semi-numerical reionisation simulation code 21cmFAST. It 
can recover constraints on astrophysical parameters from current or future 21 cm EoR experiments, accommodating 
a variety of EoR models, as well as priors on individual model parameters and the reionisation history. 
By studying the resulting impact on the EoR astrophysical constraints, 21CMMC can be used to optimise: (i) 
interferometer designs; (ii) foreground cleaning algorithms; (iii) observing strategies; (iv) alternate statistics 
characterising the 21cm signal; and (v) synergies with other observational programs. 

For more details please see http://adsabs.harvard.edu/abs/2015MNRAS.449.4246G. 

Any questions, comments or feedback (including technical support and setting up) should be directed 
to me (Bradley Greig) at bradley.greig@sns.it. If you use this code in any of your work, please cite
the above paper. 

*** 21CMMC ***

Firstly, one does not install 21CMMC on one’s computer (currently). I have modified the initial version of
CosmoHammer from the base code and so far have not changed (or tested) ‘building’
my version of CosmoHammer. Therefore clashes might occur if you already use CosmoHammer. Also, I use an
older version of CosmoHammer (the updated version can be found here https://github.com/cosmo-ethz/CosmoHammer). 
Therefore 21CMMC should run directly from the C and Python scripts included, without the need of any setup.py builds.

*** Dependencies ***

To those unfamiliar to 21cmFAST, one will need to install the FFTW3 libraries, including float support (—-enable-float).
In addition, openmp (—-enable-openmp) will also be required. Additionally, it includes standard gsl libraries.

The python code should only need scipy and numpy

*** Running 21CMMC ***

To run 21CMMC, you need to compile the C code. The folder includes a makefile for the
code. Most of the code is standard 21cmFAST C code, instead compressed into a single 
“drive_21cmMC_streamlined.c” file, which is what 21CMMC calls. It is a thinned out version 
of the EoR 21cmFAST code (mostly find_HII_bubbles.c and delta_T.c).

I have included in the folder “Boxes” z = 7,8,9 and 10 density and velocity field boxes from
21cmFAST. These were used in previous publications. Additionally, these have been used to 
generate the sensitivity curves using 21cmSense (see below). Note: These boxes were generated using WMAP 
cosmology. To use Planck cosmology, you will need to generate new boxes using 21cmFAST (the 
same is true if you want to sample the PS at a redshift other than those above).

*** Important note *** 
21CMMC only uses find_HII_bubbles.c and delta_T.c of 21cmFAST. To generate your own boxes you will need
to obtain a copy of 21cmFAST (http://homepage.sns.it/mesinger/Download.html)

The boxes included are 128^3 250 Mpc boxes. These were 768^3 boxes smoothed onto 128^3 boxes as performed
in INIT.H. These boxes are used for sampling within 21CMMC. To generate the sensitivity curves, 
I use the 21cmPS from a 1536^3 500 Mpc box, smoothed onto a 256^3 grid. Through testing, I found
these ratios to be the best. Therefore I recommend using these numbers for generating your own 
boxes for your own cosmology, redshifts or any other purposes.

*** Likelihood computation in 21CMMC ***

The likelihood computation (the main work in 21CMMC) is performed within Likelihood21cmFast.py, which is
located in “CosmoHammer_21CMMC/likelihood/module”. The likelihood computation is a chi^2 likelihood of the 21cm 
PS given some mock observation.

*** Generating Error (sensitivity curves) ***

Telescope sensitivities are not generated within 21CMMC, but rather from a separate code 
(21cmSense, https://github.com/jpober/21cmSense). The output noise files from 21cmSense are
read into 21CMMC.py and passed into the MCMC sampler. I have included some noise curves in 
“NoiseData/“. Note, these PS are binned exactly as those sampled for 250 Mpc^3 cubed boxes. At 
present there is no interpolation within 21CMMC of the 21cm PS, but this will change in a future
release.

I additionally provide the z = 8,9 and 10 sensitivity curves for a 331 tile HERA instrument with foreground wedge avoidance.

*** Important notes on 21CMMC ***

Unfortunately, 21CMMC is not completely user friendly. Owing to Python and C having to talk to one each other. I found that
the most efficient (fastest total computation time) way to run 21CMMC was to get Python to call command line code of the C 
code and then read in temporary files of the output neutral fraction and 21 cm PS.

Therefore, the current version is not very flexible and requires changes at both the Python and C interfaces to make alterations of
input variables. At the same time, there is lots of file I/O within. However, the total time of the file I/O is negligible owing to the box 
size and to the actual computation of the ionised field.

As an aside, I have developed a second version of 21CMMC (found here https://github.com/BradGreig/21CMMC_NoFileIO) 
which completely removes all file I/O within the C code, and 
removes the command line code, allowing the C code to return a value directly read into Python. However, this version is of order
50% slower than this File I/O version. This reduction in computational efficiency appears to occur owing to the way python manages its
memory. If you have experience with interfacing C and Python code or a better understanding of how to interface the 2 languages, I 
recommend looking into the other version.

Within 21CMMC.py, everything is set for running 21CMMC. Here you can set redshifts, telescope noise curves,
modelling uncertainties, parameter ranges and a few other things. Note: At present, 21CMMC is hard coded to
work only for 3 parameters (Zeta, R_mfp and Tvir) or 2 parameters (Zeta and Tvir).

*** Modifying 21CMMC ***

Unfortunately due to the hard coded nature of 21CMMC, modifying 21CMMC takes more effort that desired. However, given the base version
that I provide, it should be relatively straightforward (though slightly more effort) to implement any changes.

Basically, any changes that will need to be made must be made in drive_21cm_streamlined.c. Then, only the interface with the Python 
sampler needs to be modified to properly include the changes. Any changes to drive_21cm_streamlined.c will then need to be incorporated 
into Likelihood21cmFast.py and likely throughout 21CMMC.py, ensemble.py (CosmoHammer_21CMMC/emcee), 
CosmoHammer.py (CosmoHammer_21CMMC/sampler) and VariousInitialConditionGenerators.py (CosmoHammer_21CMMC/sampler/util)

Apologies for the completely opaque and lack of flexibility in the current version of 21CMMC!
