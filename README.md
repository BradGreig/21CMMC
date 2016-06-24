# 21CMMC
A parallelised, Monte Carlo Markov Chain analysis tool for the epoch of reionisation (EoR)

Latest release: 24/06/2016
- Added support of observational priors
- Added a fourth parameter (power-law index) for mass-dependent ionising efficiency
- Mock observation and noise files no longer need to be binned the same (now interpolated)

21CMMC is a fast and efficient Python sampler of the semi-numerical reionisation simulation code 21cmFAST. It 
can recover constraints on astrophysical parameters from current or future 21 cm EoR experiments, accommodating 
a variety of EoR models, as well as priors on individual model parameters and the reionisation history. 
By studying the resulting impact on the EoR astrophysical constraints, 21CMMC can be used to optimise: (i) 
interferometer designs; (ii) foreground cleaning algorithms; (iii) observing strategies; (iv) alternate statistics 
characterising the 21cm signal; and (v) synergies with other observational programs. 

For more details please see http://adsabs.harvard.edu/abs/2015MNRAS.449.4246G. For a more detailed readme, open readme.txt
within the 21CMMC_SourceCode directory (in Programs)

Any questions, comments or feedback (including technical support and setting up) should be directed 
to me (Bradley Greig) at bradley.greig@sns.it. If you use this code in any of your work, please cite
the above paper and link to this repository. 
