---
title: 'Electronic Supplementary Materials: Analyses workflow of Vea and Grimaldi feeding Coccoidea (Hemiptera) in the light of amber inclusions' : A time scale for scales: divergence time and diversification rates of the plant-sap
output: html_document
---


#Summary
This supplementary material document details the analyses performed for the results presented in Vea and Grimaldi (submitted). A Time scale for scales: divergence time and diversification rates of the plant-sap feeding Coccoidea (Hemiptera) in the light of amber inclusions. 
It includes command lines used in the MrBayes vers. 3.2.3 analyses, as well as command lines for the R analyses. Finally, this file also includes supplementary figures referred in the main paper.

#Data

  - Combined dataset downloadable in Dryad (url to be added after manuscrip accepted): This dataset includes XX taxa, with 174 morphological characters and xx aligned nucleotides including 18S, 28S and EF-1, in nexus format.
  
The dataset includes the following command line block, defining the analysis basic settings:

```
begin mrbayes;
	set autoclose=yes nowarn=yes;
	[Characters sets were defined with male and female characters as separate sets.]
	[A subset for EF1 is defined with 1st and 2nd codon position together and separated from 3rd codon]
	charset MorphoM = 1-126;
	charset MorphoF = 127-174;
	charset 18S = 175-882;
	charset 28SD2D3 = 883-1846;
	charset 28SD10 = 1847-2903;
	charset EF-1a = 2904-3905;
	charset EF-1a_12 = 2904-3905\3 2905-3905\3;
	charset EF-1a_1 = 2904-3905\3;
	charset EF-1a_2 = 2905-3905\3;
	charset EF-1a_3 = 2906-3905\3;
	[From these character sets, 6 different partitions were separated]
	partition seven = 7: MorphoM, MorphoF, 18S, 28SD2D3, 28SD10, EF-1a_12, EF-1a_3;
	partition eight = 8: MorphoM, MorphoF, 18S, 28SD2D3, 28SD10, EF-1a_1, EF-1a_2, EF-1a_3;
	partition six = 6: MorphoM MorphoF, 18S, 28SD2D3, 28SD10, EF-1a_12, EF-1a_3;
	partition five = 5: MorphoM MorphoF, 18S, 28SD2D3 28SD10, EF-1a_12, EF-1a_3;
	
	[Five partitions were used in this analysis]
	set partition = five;	
	
	[Set evolutionary models based on ModelTest]

	lset applyto=(1)   coding=variable rates=gamma; [Mk model for morphology]
	lset applyto = (2,3,4) nucmodel=4by4 nst=6 rates=invgamma covarion=no;[GTR+I+G for 18S, 28S and EF1a1]
	lset applyto=(5) nst=2 rates=invgamma; [HKY+G+I for EF1a 2 and 3]
	
	unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all); 
	prset applyto=(all) ratepr=variable;
	
	showmodel; [showing models that we just set up]
	
	[Outgroup was defined as an Aphididae]
	outgroup Acyrthosiphon_pisum;
	
	[Setting fossil set]

	taxset fossils= Alacrena_peculiaris   	
Albicoccus_dimai   			
Apticoccus_fortis   		
Apticoccus_longitenuis   	
Apticoccus_minutus   		
undescribed_ARC60_1   					
Arnoldus_capitatus   		
Burmacoccus_danyi   		
Cretorthezia_hammanaica   	
Electrococcus_canadensis   	
Eomatsucoccus_casei   		
Grimaldiella_gregaria   	
Grohnus_eichmanni   		
Heteromargarodes_hukamsinghi 
Hodgsonicoccus_patefactus   
Inka_minuta   				
Jersicoccus_kurthi   		
Gilderius_eukrinops   		
Kozarius_achronus   	
Kozarius_perpetuus   	
Kuenowicoccus_pietrzeniukae  
Kukaspis_usingeri   		
Lebanococcus_longiventris   
Magnilens_glaesaria   	
Marmyan_barbarae   			
Normarkicoccus_cambayae   	
Palaeosteingelia_acrai   	
Palaeotupo_danieleae   		
Palaeonewsteadia_huaniae   	
Pedicellicoccus_marginatus   
Pennygullania_electrina   	
Pityococcus_moniliformalis   	
Protorthezia_aurea   		
Pseudoweitschatus_audebertis 
Rosahendersonia_prisca   	
Serafinus_acutipterus   	
Solicoccus_nascimbenei   	
Steingelia_cretacea   		
Turonicoccus_beardsleyi   	
Weitschatus_stigmatus   	
Williamsicoccus_megalops   	
Xiphos_vani   				
Xylococcus_granbenhorstii; 

	[Defining family constraints]
constraint root = 1-.;
constraint coccoidea_nofossils -1 = 6-78;
constraint coccoidea_wfossils -1 = 6-121;
constraint coccidae -1 = 13 17 32 113;
constraint diaspididae -1 = 8 15 56 104;
constraint kermesidae -1= 41 42;
constraint kuwaniidae -1= 43 51;
constraint margarodidae -1 = 24 25 33 37 62 92;
constraint matsucoccidae -1 = 46 47 48 49 50 89;
constraint ortheziidae -1 = 40 52 53 63 87 107;
constraint pseudococcidae -1= 7 14 16 29 34 61 64 65 96 119;
constraint putoidae -1 = 66 67 68 69 70;
constraint xylococcidae -1 = 77 78 121;

```

#Analyses
All MrBayes 3.2.X analyses were perfomed through the Cipres Gateway Portal (url) with four replicate runs of 10-20 million generations.

##Preliminary phylogenetic analyses 
Using the combined dataset nexus file, the following were added command lines within the `begin mrbayes;` block (respectively for each analysis).

1. Non-clock analysis

```
delete fossils;
mcmcp temp=0.2 nchain=4 samplefreq=1000 printfr=100 nruns=4;
mcmcp filename=nonclock;
mcmc ngen=20000000;
```

2. Strict-clock analysis

```
delete fossils;
prset brlenspr=clock:uniform;
mcmcp temp=0.2 nchain=4 samplefreq=1000 printfr=100 nruns=4;
mcmcp filename=strictclock;
mcmc ngen=10000000;
```
After the analyses were finished, output files were downloaded and loaded into local MrBayes to check for convergence (avergae standard deviation of split frequency < 0.05), and summary topology base on the all compatibility trees were generated as follow:

```
sump nruns=4;
sumt filename=nonclock nruns=4 contype=allcompat;
```

After performing `sumt` for both preliminary analyses, the following outputs were used for next analyses: 
  
  - nonclock.con.tre
  - strictclock.con.tre

##IGR prior values
To obtain the value for the IGR model, we used the R script from Ronquist et al. (2012). The original script had to have the outgroup replaced to Acyrthosiphon_pisum. For details on the analysis performed in R deposited in Github [here] (url to github IGR). The prior values obtained to set the IGR model are summarized in Table SX of Electronic Supplementary Material.


##Calibrated analyses
The following command lines were added to set the IGR model priors for both node-dating and tip-dating analyses.


```
[values obtained with unconstrained nonclock and strict clock analyses]
delete fossils;
prset brlenspr=clock:uniform;
prset clockvarpr=igr;	
prset igrvarpr=exp(25.64493156); [ 1/x with x=median variance obtained from R script analysis comparing non clock and strict clock branch lengths]
prset nodeagepr = calibrated;
prset clockratepr = lognorm(-6.060527447,0.051884487); 

```


Eleven different calibrated analyses were run in MrBayes with different calibration prior combinaisons. Information on calibrations priors are summarized in Table S3.


1. ND-A: all 12 calibrations with the following command line in MrBayes 3.2.3


```
prset topologypr=constraints(root,coccoidea_nofossils,coccidae,diaspididae, kermesidae, kuwaniidae,margarodidae,matsucoccidae,ortheziidae,pseudococcidae,putoidae,xylococcidae);

[this line was modified according to analysis ND-A, ND-Anoroot or ND-Alognormal]
calibrate root=offsetexp(240,250); 
[calibrate root=lognormal(240,10); ]


calibrate coccoidea_nofossils=offsetexp(140, 240);[2-definitive coccoidea min: 140 not likely to be more than 240]
calibrate coccidae=offsetexp(98,110); [3]
calibrate diaspididae=offsetexp(50,100);[4]
calibrate kermesidae=offsetexp(45,100); [5]
calibrate kuwaniidae=offsetexp(45,100); [6]
calibrate margarodidae=offsetexp(50,100); [7]
calibrate matsucoccidae=offsetexp(135,140);	[8]	
calibrate ortheziidae=offsetexp(135,140);	[9]	
calibrate pseudococcidae=offsetexp(135,140); [10]
calibrate putoidae=offsetexp(45,100); [11]
calibrate xylococcidae=offsetexp(135,140); [12]
```

2. ND-B: excluding calibration 3
3. ND-C: excluding calibration 4
4. ND-D: excluding calibration 5
5. ND-E: excluding calibration 6
6. ND-F: excluding calibration 7
7. ND-G: excluding calibration 8
8. ND-H: excluding calibration 9
9. ND-I: excluding calibration 10
10. ND-J: excluding calibration 11
11. ND-K: excluding calibration 12

##Total evidence analyses
Four total evidence analyses were performed with 43 prior ages set as fixed for fossil terminals (see Table S2 for fossil ages), as well as the alternative options:

1. TD-A: 43 fixed age priors at fossil tips, and a node prior set at the root following an off set exponential distribution
2. TD-B: 43 fixed age priors at fossil tips only
3. TD-C: 43 fixed age priors at fossil tips, and a node prior set at the root following a lognormal distribution


The following command lines were added, replacing the node calibration settings from above, with last two command lines removed according to analysis:

```
calibrate   Albicoccus_dimai=fixed(98)
  			Apticoccus_minutus=fixed(135) 
			Apticoccus_fortis=fixed(135) 	
			Arnoldus_capitatus=fixed(45) 
			Burmacoccus_danyi=fixed(98) 
			Cretorthezia_hammanaica=fixed(135) 
			Eomatsucoccus_casei=fixed(92) 
			Electrococcus_canadensis=fixed(100) 
			Grimaldiellia_gragaria=fixed(92) 
			Grohnus_eichmanni=fixed(45) 
			Inka_minuta=fixed(85) 
			Jersicoccus_kurthi=fixed(92)
			Kuenowicoccus_pietzeniukae=fixed(45) 
			Kukaspis_usingeri=fixed(100) 
			Lebanococcus_longiventris=fixed(135) 
			Hodgsonicoccus_patefactus=fixed(135) 
			Xiphos_vani=fixed(135) 
			Apticoccus_longitenuis=fixed(135) 
			Marmyan_barbarae=fixed(98)
			Palaeosteingelia_acrai=fixed(135) 
			Palaeotupo_danielae=fixed(135) 
			Palaeonewsteadia_huaniae=fixed(45) 
			Pennygullania_electrina=fixed(135) 
			Pityococcus_moniliformalis=fixed(45) 
			Protorthezia_aurea=fixed(45) 
			Serafinus_acupiterus=fixed(45)
			Solicoccus_nascimbenei=fixed(92) 
			Steingelia_cretacea=fixed(92)
			Heteromargarodes_hukamsinghi=fixed(50) 
			Normarkicoccus_cambayae=fixed(50)
			Kozarius_perpetuus=fixed(98) 
			Kozarius_achronus=fixed(98) 
			Gilderius_eukrinops=fixed(98) 
			Pseudoweitschatus_audebertis=fixed(98) 
			Rosahendersonia_prisca=fixed(98) 
			Magnilens_glaesaria=fixed(98) 
			Pedicellicoccus_marginatus=fixed(98) 
			Alacrena_peculiaris=fixed(98) 
			Williamsicoccus_megalops=fixed(135)
			undescribed_ARC60_1=fixed(100) 
			Turonicoccus_beardlseyi=fixed(92)
			Weitschatus_stigmatus=fixed(45) 
			Xylococcus_grabenhorstii=fixed(45)
			;   


[Presettings for the total evidence analysis]

prset brlenspr=clock:uniform;
prset clockvarpr=igr;	
prset igrvarpr=exp(25.64493156); [ 1/x with x=median variance obtained from R script analysis comparing non clock and strict clock branch lengths]
prset nodeagepr = calibrated;
prset clockratepr = lognorm(-6.060527447,0.051884487); 

prset topologypr=constraints(root,coccoidea_wfossils);

[removed according to one of three types of analyses TD-A, TD-B or TD-C]
calibrate root=offsetexp(240,250);
[calibrate root=lognormal(240,10);]

```
MCMC for Node-dating and Tip-dating analyses

```
mcmcp temp=0.1 nchain=4 samplefreq=1000 printfr=100 nruns=4;
mcmcp filename=filename;
mcmc ngen=40000000;
```


Tree convergence were checked and summary tree were compiled using the all compatibility summary option.

```
sump nruns=4;
sumt contype=allcompat;
end;
```



