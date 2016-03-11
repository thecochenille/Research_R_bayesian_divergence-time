------------------------------
# Analyses worflow for the manuscript Vea and Grimaldi: 
##Putting scales into evolutionary time: the divergence of major scale insect lineages (Hemiptera) predates the radiation of modern angiosperm hosts


####Scientific Reports 6, 23487 doi: 10.1038/srep23487
####Repository status: public

-----------------------------



# Summary
This repository includes the comman lines used for MrBayes analyses (divergence time estimates) and other analyses performed in R (IGR model estimation and LTT plots). Only analyses from R have a direct output using Rmd files.


#Files included
1. Nexus files including all topologies obtained with the MrBayes command lines presented in this repository. These are:
	- `nonclock.con.tre`: from preliminary MrBayes analysis to be used in `Coccomorpha-IGR.Rmd`
	- `strictclock.con.tre`: from preliminary MrBayes analysis to be used in `Coccomorpha-IGR.Rmd`
	- `ND-A-offsetexp.tre` : used in `LTTplots.Rmd`
	- `ND-A-lognormal.tre`: used in `LTTplots.Rmd`
	- `ND-A-noroot.tre`: used in `LTTplots.Rmd`
	- `ND-B.tre`: used in `LTTplots.Rmd`
	- `ND-B.tre`: used in `LTTplots.Rmd`
	- `ND-B.tre`: used in `LTTplots.Rmd`
	- `ND-B.tre`: used in `LTTplots.Rmd`
	- `ND-B.tre`: used in `LTTplots.Rmd`
	- `ND-B.tre`: used in `LTTplots.Rmd`
	- `ND-B.tre`: used in `LTTplots.Rmd`
	- `ND-B.tre`: used in `LTTplots.Rmd`
	- `ND-B.tre`: used in `LTTplots.Rmd`
	- `ND-B.tre`: used in `LTTplots.Rmd`
	- `TD-offsetexp.tre`: used in `LTTplots.Rmd`
	- `TD-lognormal.tre`: used in `LTTplots.Rmd`
	- `TD-noroot.tre`: used in `LTTplots.Rmd`


2. `Coccomorpha-IGR.md`: workflow and modified script to estimate IGR model (the necessary topologies from preliminary analyses are described in the `divergence-time.Rmd` file below).

3. `divergence-time.md`: workflow of MrBayes analyses and LTT plots in R.

4. `LTTplots.md`: workflow and R script to obtain the LTT plots discussed in the Electronic Supplementary Material. 

5. Folder `Datasets` including all datasets from unaligned sequences to combined dataset ready for analysis


#Where to start?
Although each file can be viewed individually by following the manuscript and EMS, the best way to examine the worklow of the analyses is to start reading [divergence-time.md](https://github.com/zourloubidou/Coccomorpha-divergence-time/blob/master/divergence-time.md).

#Reference
Vea I. and D. Grimaldi. (2016) Putting scales into evolutionary time: the divergence of major scale insect lineages (Hemiptera) predates the radiation of modern angiosperm hosts. Scientific reports 6, 23487 doi: 10.1038/srep23487
