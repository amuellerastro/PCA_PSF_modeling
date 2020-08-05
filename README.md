The routine to be called is PCA_LOCI_ADI.pro.

It depends on several subroutines provided by different IDL libraries: astrolib, coyote, mpfit, StarFinder, textoidl, mperrin.

It only works for data products produced by my pipeline (batch_NACO_reduce.pro)!


-data products (image registered cubes, PSF, and parallactic angles) have to be following the definition (naming, structure, location) set by my NACO pipeline (batch_NACO_reduce.pro)

-line 341: change locatio of the directory of datapaths.txt. This text file contains the absolute path to the reduced data. One entry looks like e.g. /home/amueller/NACO/data/HD123/

-The data to be processed can only be fro the spectral range Lp, Mp or Ks.

-The strength of this routine is the - as I call it - "dynamic smart PCA". In the terminal you will be asked to select the method to be applied. For this method select "v".

-The routines allow to perform a recentering of the image cube using amoeba. I am not recommending to use it. Select "n" in the terminal.

-The LLSG method is not fully implemented and should not be applied.

-line 933: modify the path pointing to the transmission curves of the Ks and ND_short filter.

-line 535: candidate detection: if set to '1' automatic candidate detection is applied. It is experimental and not robust close to the IWA. 

-line 536: fake planet injection: if set to '1' fake planets are injected and contrast curves are produced. But increases computational time.

-line 537: apply spatial filtering if set to '1'. There are 9 different filter options which can be selected. It is not recommended to use this feature since it changes the contrast significantly and is not well accepted in the community.

-Line 540-548: if you are using "v" (the preferred method) you can define a range of values for each parameter. The current values are the most useful ones. It will produce dozens of PCA reduced images and allow to identify the best possible parameter setup. The individual parameters are defined in the name of the result directory (e.g. HD123/PCA_Reduction/sPCAvar_noFilter_Step1.0_Rin0.5_Delta0.5/) and in setup_Lp.txt located inside the individual directories.
  -section: experimental, set to 0
  -rin_init: the size of the central mask applied to the image in fraction of FWHM
  -klip: number of pca modes to be used
  -delta: protection angle in fractions of FWHM
  -step_init: width of the annuli in FWHM. No need to change this.

-line 619: allows you to deinfe a new image size by cutting the image around the image center.

-There are dozens of intermediate products written out. Their file name explains the meaning. The most important files are:
  -img_Lp_rsz_scaled_pca_median.fits  The final PCA reduced image.
  -plot_Lp_contrast_map.pdf 5-sigma contrats map (if fcpflag='1')
  -Lp_SNR_map.fits  SNR map (if fcpflag='1')
  -<id>_Lp_<night>_overview_PCAvar.pdf  Summary of the observation (if fcpflag='1')
  -plot_Lp_contrastMag.ps plot of 5-sigma contrast (if fcpflag='1')
  
-no correction for true north is applied. AGPM throughput correction is applied.
