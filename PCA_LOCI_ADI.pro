@merging_window_v5.pro
@readfits.pro
@sxpar.pro
@gettok.pro
@valid_num.pro
@writefits.pro
@check_fits.pro
@fxpar.pro
@fxaddpar.pro
@sxdelpar.pro
;@pca_adi_v38.pro
@pca_adi_v39.pro
@pca_adi_v39_am.pro
@dist.pro
@mean.pro
@moment.pro
@reverse.pro
@rot.pro
@mkhdr.pro
@sxaddpar.pro
@get_date.pro
@daycnv.pro
; @contrast_curve_v71.pro
@filter_image.pro
@factor.pro
@prime.pro
@psf_gaussian.pro
@gaussian.pro
@convolve.pro
; @convol_fft.pro
@real_part.pro
@cv_coord.pro
@robust_sigma.pro
@bandpass_filter_MOD.pro	;modification w.r.t. to FFT because IDL<7.1 had no /center keyword
@mfilter_v4b.pro
@spca_adi_v13.pro
@spca_adi_v15.pro
@congrid.pro
; @shifti.pro
; @fcp_v83.pro
; @alpha_r_v61.pro
@interpol.pro
@image_stddev.pro
; @prep_ps.pro
@loadct.pro
; @filepath.pro
; @path_sep.pro
; @tvim.pro
; @color_key.pro
; @end_ps.pro
@sigfig.pro
@mrdfits.pro
@fxposit.pro
@mrd_hread.pro
@mrd_skip.pro
@cgcolor.pro
@cggetcolorstate.pro
@strsplit.pro
@cgimage.pro
@image_dimensions.pro
@setdefaultvalue.pro
@cgdefaultcolor.pro
@ps_background.pro
@cgcolorfill.pro
@cgsetcolorstate.pro
@cgcheckforsymbols.pro
@cgplot.pro
@cgbitget.pro
@convert_to_type.pro
@colorsareidentical.pro
@cgcolorbar.pro
@cgdefcharsize.pro
@array_indices.pro
@fcp_am.pro
@noise_curve_am.pro
@contrast_curve_am.pro
@sph_ird_transmission.pro
@strnumber.pro
@circint_MOD.pro
@proceeding_text.pro
@sumann_MOD.pro
@readcol.pro
@remchar.pro
@linspace.pro
@cgetrng.pro
@pixwt.pro
@ts_diff.pro
@mfilter_am.pro
@logm.pro
@headfits.pro
@get_eso_keyword.pro
@shift_sub.pro
@alpha_r_am.pro
@robust_mean.pro
@avg.pro
@loci_adi_v19.pro
@nnls.pro
@bvls.pro
@adi_v01.pro
@transmatch9.pro
@whereprams.pro
@detection_limit_small_sample.pro
@isarray.pro
@gauss_pdf.pro
@t_cvf.pro
@t_pdf.pro
@ibeta.pro
@bisect_pdf.pro
@aper.pro
@legend.pro
@cleanplot.pro
@strn.pro
@textoidl.pro
@textable.pro
@strtrans.pro
@translate_sub_super.pro
@nexttok.pro
@closest2.pro
; @planet_lower_mass_limit.pro
@min_curve_surf.pro
@cgsurf.pro
@cgsnapshot.pro
@cgcolor24.pro
@t3d.pro
@cgplots.pro
@cgsymcat.pro
@dist_circle.pro
@explabel.pro
@setdecomposedstate.pro
@decomposedcolor.pro
@get_w.pro
@amoeba.pro
@display.pro
@imgexp.pro
@imgscl.pro
@rmsmap.pro
@rmsmap_am.pro
@fftshift.pro
@boot_mean.pro
@llsg_am.pro
@define_annuli.pro
@get_annulus_quad.pro
@patch_rlrps.pro
@thresholding.pro
@meanabsdev.pro
@sign.pro
@caldat.pro
@roundness.pro
@psf.pro
@nint.pro
@rangegen.pro
@candidate_detection.pro
@bsort.pro
@spca_adi_var.pro
@spca_adi_var_corr.pro
@clipscl.pro
@cgscalevector.pro
@fpufix.pro
@showsym.pro
@ten.pro
@cgaxis.pro
@cgplot.pro
@cgoplot.pro
@parangle.pro
@gauss_noise_std.pro
@median_filter.pro
@size52.pro
@histo.pro
@histo_hwhm.pro
@gaussfit.pro
@poly_fit.pro
@curvefit.pro
@starfinder.pro
@fwhm.pro
@get_max.pro
@subs_to_coord.pro
@peak_width.pro
@peak_area.pro
@binary_array.pro
@image_core.pro
@search2d.pro
@sub_array.pro
@array_overlap.pro
@shifted_templates.pro
@image_shift.pro
@extend_array.pro
@search_objects.pro
@all_max.pro
@starlist.pro
@star.pro
@update_list.pro
@sort_list.pro
@merge_list.pro
@star_param.pro
@max_search.pro
@distance.pro
@correlate_max.pro
@extract_overlap.pro
@correlation_coeff.pro
@where_stars.pro
@add_subscript.pro
@fitstars.pro
@frequency.pro
@newton_gauss.pro
@image_model.pro
@add_overlap.pro
@plane.pro
@ls_sys.pro
@diag_mult.pro
@min_norm_inversion.pro
@ginv.pro
@convergence.pro
@relative_error.pro
@fitting_errors.pro
@compare_lists.pro
@reciprocal_distance.pro
@extract_stars.pro
@cgerase.pro
@cgresizeimage.pro
@closest.pro
@sky.pro
@mmm.pro
@asinh.pro
@scale_image_am.pro
@fits_add_checksum.pro
@checksum32.pro
@n_bytes.pro
@is_ieee_big.pro
@host_to_ieee.pro
@fits_ascii_encode.pro
@correl_images.pro
@fxmove.pro
@wiener.pro
@sadi.pro

pro PCA_LOCI_ADI, sphere=sphere, naco=naco, nici=nici, lmir=lmir, magao=magao, nirc2=nirc2, right_handed=right_handed, pmodel=pmodel, plotonly=plotonly, cut=cut, disk=disk
;cut: cut out bright CC for contrast computation
;disk: for sPCAvar fix KLIP to 5

;pmodel: if set then btsettl models used to estimate planet mass - NOW called planet_mass_detection_limit.pro and is a stand alone routine

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

if keyword_set(plotonly) then begin

  print, ''
  print, 'Parameters have to match!'
  print, ''

endif

if keyword_set(sphere) then selins = 's'
if keyword_set(naco) then selins = 'n'
if keyword_set(magao) then selins = 'm'
if keyword_set(nirc2) then selins = 'n2'

if not keyword_set(naco) and not keyword_set(sphere) and not keyword_set(nici) and not keyword_set(lmir) and not keyword_set(magao) then begin

  selins = ''
  read, 'Select Instrument Naco (n) / Sphere (s) / NICI (ni) / LMIRCam (l) / MagAo (m) / Nirc2 (n2): ', selins
; selins = 'n'
  if (selins eq 'n') then naco = 1
  if (selins eq 'ni') then nici = 1
  if (selins eq 's') then sphere = 1
  if (selins eq 'l') then lmir = 1
  if (selins eq 'm') then magao = 1
  if (selins eq 'n2') then nirc2 = 1

endif

;========================================================================

;select star and get its needed properties

filestar = file_search('/home/amueller/work/IDLlibs/AO/TargetProperties/Targets/*.sav', count=nstars)
stars = strarr(nstars)
for i=0,nstars-1 do begin

  stars[i] = strmid(filestar[i], 56, strlen(filestar[i])-56)
  stars[i] = strmid(stars[i], 0, strlen(stars[i])-4)

endfor

nstars = n_elements(stars)
print, ''
for i=0,nstars-1 do print, strcompress(i+1, /rem), ' ', stars[i]
read, 'Select Star: ', selstar
; selstar = '109'
star = stars[selstar-1]
restore, filestar[selstar-1], /verbose

if keyword_set(sphere) then stellmag = [jmag, hmag, kmag]	;JHK
if keyword_set(naco) then stellmag = [lmag, kmag]	;L'
if keyword_set(nici) then stellmag = [lmag]	;L'
if keyword_set(lmir) then stellmag = [lmag,lmag]	;L'

;========================================================================

if keyword_set(sphere) then begin

  readcol, '/home/amueller/work/IDLlibs/AO/SPHERE/datapaths.txt', tmp, format='a', /silent
  print, ''

  idx = where(strmatch(tmp, '*'+star+'*') eq 1)
  tmp = tmp[idx]

  for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]
  print, ''
  read, 'Select path: ', selp
  datadir = tmp[selp-1]+'/IRDIS/'	;reduced data
  path = tmp[selp-1]+'/RAW/IRDIS/'	;raw data

  if (strmatch(datadir, '*'+star+'*') ne 1) then begin

    print, ''
    print, 'Wrong star selected! Stop.'
    stop

  endif

  ;selection of used filter pair
  print, ''
  filter_pair = ['K1/K2', 'H2/H3', 'J2/J3', 'BB_Y', 'BB_J', 'BB_H', 'BB_K', 'PaB/PaB']
  for i=0,n_elements(filter_pair)-1 do print, strcompress(i+1,/rem), ' ', filter_pair[i]
  read, 'Select filter pair: ', self
  if (self eq 1) then camera_filter = ['K1', 'K2']
  if (self eq 2) then camera_filter = ['H2', 'H3']
  if (self eq 3) then camera_filter = ['J2', 'J3']
  if (self eq 4) then camera_filter = ['YL', 'YR']
  if (self eq 5) then camera_filter = ['JL', 'JR']
  if (self eq 6) then camera_filter = ['HL', 'HR']
  if (self eq 7) then camera_filter = ['KL', 'KR']
  if (self eq 8) then camera_filter = ['PaB']
  
endif

;========================================================================

if keyword_set(naco) then begin

  readcol, '/home/amueller/work/IDLlibs/AO/NACO/datapaths.txt', tmp, format='a', /silent
  print, ''

  idx = where(strmatch(tmp, '*'+star+'*') eq 1)
  tmp = tmp[idx]

  for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]
  print, ''

  read, 'Select path: ', selp
; selp = 1
  datadir = tmp[selp-1]	;reduced data
  path = tmp[selp-1]+'/RAW/'	;raw data

  if (strmatch(datadir, '*'+star+'*') ne 1) then begin

    print, ''
    print, 'Wrong star selected! Stop.'
    stop

  endif
  
  ;selection of used filter
  print, ''
  filter_pair = ['Lp', 'Mp', 'Ks']
  for i=0,n_elements(filter_pair)-1 do print, strcompress(i+1,/rem), ' ', filter_pair[i]
  read, 'Select filter pair: ', self
  if (self eq 1) then camera_filter = ['Lp']
  if (self eq 2) then camera_filter = ['Mp']
  if (self eq 3) then camera_filter = ['Ks']

endif

;========================================================================

if keyword_set(nici) then begin

  readcol, '/home/amueller/work/IDLlibs/AO/NICI/datapaths.txt', tmp, format='a', /silent
  print, ''

  idx = where(strmatch(tmp, '*'+star+'*') eq 1)
  tmp = tmp[idx]

  for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]
  print, ''

  read, 'Select path: ', selp
; selp = 1
  datadir = tmp[selp-1]	;reduced data
  path = tmp[selp-1]+'/RAW/'	;raw data

  if (strmatch(datadir, '*'+star+'*') ne 1) then begin

    print, ''
    print, 'Wrong star selected! Stop.'
    stop

  endif

  camera_filter = ['Lp']

endif

;========================================================================

if keyword_set(nirc2) then begin

  readcol, '/home/amueller/work/IDLlibs/AO/NIRC2/datapaths.txt', tmp, format='a', /silent
  print, ''

  idx = where(strmatch(tmp, '*'+star+'*') eq 1)
  tmp = tmp[idx]

  for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]
  print, ''

  read, 'Select path: ', selp
; selp = 1
  datadir = tmp[selp-1]	;reduced data
  path = tmp[selp-1]+'/RAW/'	;raw data

  if (strmatch(datadir, '*'+star+'*') ne 1) then begin

    print, ''
    print, 'Wrong star selected! Stop.'
    stop

  endif
  
  ;selection of used filter
  print, ''
  filter_pair = ['Lp', 'Ms']
  for i=0,n_elements(filter_pair)-1 do print, strcompress(i+1,/rem), ' ', filter_pair[i]
  read, 'Select filter pair: ', self
  if (self eq 1) then camera_filter = ['Lp']
  if (self eq 2) then camera_filter = ['Ms']

endif

;========================================================================

if keyword_set(lmir) then begin

  readcol, '/home/amueller/work/IDLlibs/AO/LMIRCam/datapaths.txt', tmp, format='a', /silent
  print, ''

  idx = where(strmatch(tmp, '*'+star+'*') eq 1)
  tmp = tmp[idx]

  for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]
  print, ''
  read, 'Select path: ', selp
  datadir = tmp[selp-1]+'/'	;reduced data
  path = tmp[selp-1]+'/RAW/'	;raw data

  if (strmatch(datadir, '*'+star+'*') ne 1) then begin

    print, ''
    print, 'Wrong star selected! Stop.'
    stop

  endif

  ;selection of used filter pair
  camera_filter = ['dx', 'sx']

endif

if keyword_set(magao) then begin

  readcol, '/home/amueller/work/IDLlibs/AO/MagAO/datapaths.txt', tmp, format='a', /silent
  print, ''

  idx = where(strmatch(tmp, '*'+star+'*') eq 1)
  tmp = tmp[idx]

  for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]
  print, ''
  read, 'Select path: ', selp
  datadir = tmp[selp-1]+'/'	;reduced data
  path = tmp[selp-1]+'/RAW/'	;raw data

  if (strmatch(datadir, '*'+star+'*') ne 1) then begin

    print, ''
    print, 'Wrong star selected! Stop.'
    stop

  endif

  ;selection of used filter pair
  camera_filter = ['Ha', 'CntHa']

endif

;========================================================================

; methods = ''
; read, 'variable sPCAcorr / sPCAvar / sPCA / PCA Box / LLSG / LOCI / ADI / smartADI (vc/v/p/c/l/loci/a/sa): ', methods
methods = ['v']
; methods = ['p']
; methods = ['p'];,'l']

; recenter = ''
; read, 'Recenter images (y/n): ', recenter
recenter = 'n'
;if (recenter eq 'n') then rcflag = ''

for mth=0,n_elements(methods)-1 do begin

  method = methods[mth]

  if (method eq 'c') then dirflag = '_PCAbox'
  if (method eq 'vc') then dirflag = '_PCAvarcorr'
  if (method eq 'v') then dirflag = '_PCAvar'
  if (method eq 'p') then dirflag = '_PCA'
  if (method eq 'l') then dirflag = '_LLSG'
  if (method eq 'a') then dirflag = '_ADI'
  if (method eq 'sa') then dirflag = '_sADI'

  if (method eq 'v' or method eq 'vc') then begin

    ;pointdisk = ''
    ;read, 'Point source or disk (p/d): ', pointdisk
    if keyword_set(disk) then pointdisk = 'd' else pointdisk = 'p'

  endif else begin

    pointdisk = ''

  endelse

  ;========================================================================

  canddet = '0'	;automatic candidate detection, 0=off
  fcpflag = '0'	;with or without fake companion injection
  filterflag = ['0'];['0', '1']	;spatial filtering applied before PCA reduction
				    ;if lowpass filter is changed, adjust gaussian smoothing in mfilter_4b and mfilter_am
				    ;if filter method is changed, don't forget to change mfilter_am keyword in fcp_am for PSF filtering
  if (method eq 'v' or method eq 'vc') then begin

    section = 0	;1: divide image in radial anuli with sections / 0: divide image in annuli only
    rin_init = [0.1,0.25,0.5,0.75,1.];[0.5]
    klip = uint([5,10,20,30]);uint([20])
    delta = [0.1,0.2,0.3,0.5,0.75,1.];[0.75]
    step_init = [1]

  endif

  if (method eq 'c') then begin

    zone = [11];, 100, 50, 25, 15]	;square zone on which the PCA will be perfomed [in pixels]
					  ;note: an odd number of boxes should better be used if star centred in the images)
    ;maximum value of klip is zone^2, e.g. zone = 5 -> max klip=25
    klip = uint([50])
    delta = [1.0]
    rin_init = [1]

  endif

  if (method eq 'p') then begin

    rin_init = [2.];,3.];, 4]	;inner radius where PCA starts, in resolution elements (fwhm) 
    klip = uint([5,10,15,20]);,15,20,25];,20]	;number of principal components to be used in the (s)PCA
    step_init = [1];, 1]	;width of the annuli for sPCA, in resolution elements (fwhm)
    delta = [0.25, 0.3, 0.35, 0.4, 0.5];[0.25,0.5];,1.0]	;in FWHM, parallactic angle exclusion region for building the PSF library in sPCA mode, defined as the minimum angular separation in units of $\lambda/D$ between current companion position other images of the cube. Default = 1*lambda/D.

  endif

  ;if (method eq 'p') then smethod = '3';read, 'v13 or v15 (3/5): ', smethod
  smethod = '3'


  if (method eq 'a' or method eq 'c') then rin_init = [0.1,0.5,1.]

  if (method eq 'sa') then begin

    rin_init = [0.1,0.25,0.5,1.]
    delta = [0.5,1.0,2.0]
    distdelta = 0.2	;distance in arcsec where protection is computed

  endif

  if (method eq 'l') then begin

    klip = uint([10])	;originally called rank
    thresh = 1.
    max_iter = 10.
    low_rank_mode = 'svd'
    thresh_mode = 'soft'
    nzone = 1
    delta = [0]	;so far no frame selection implemented. Not possible?
    rin_init = [0]

  endif


  if (method eq 'loci') then begin

    rin_init = [1.];,3.];, 4]	;inner radius where PCA starts, in resolution elements (fwhm) 
    klip = [16]
    step_init = [1]
    delta = [0.4]
    Na = 250.	;number of resolution element inside the optimization zone, too small can kill companion signal
    g = 2.0	; geometry of the optimization zone (see Lafreniere et al. 2007)
		  ; 0.75* = elongated along theta, squashed radially
		  ; 1     = sqr
		  ; 2     = elongated radially, 2xlong
  endif


  if keyword_set(right_handed) then right_handed = 1 else right_handed = 0
  if keyword_set(nici) then right_handed = 1
  if keyword_set(magao) then right_handed = 1
  ;read, 'Right Handed 0 (science) / 1 (blind test): ', right_handed

  hdr = headfits(datadir+'img_'+camera_filter[0]+'_dc.fits', exten=0, /silent)
  dim = float(get_eso_keyword(hdr, 'NAXIS1'))
  ;dim = 300.	;cuting image down to this size, EVEN number

  ;check if the reduced cubes have header info as implemented later in the reduction process
  ;if file is older then we have to select one science and one flux frame
  if (selins eq 'n') then tmp = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DET DIT'),/rem)
  if (selins eq 's') then tmp = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO DET SEQ1 DIT'),/rem)
  if (selins eq 'ni') then tmp = ''
  if (selins eq 'l') then tmp = strcompress(get_eso_keyword(hdr, 'ITIME'),/rem)
  if (selins eq 'm') then tmp = strcompress(get_eso_keyword(hdr, 'DIT'),/rem)
  if (selins eq 'n2') then tmp = strcompress(get_eso_keyword(hdr, 'ITIME'),/rem)
  if (tmp eq '') then flagnewred = 0 else flagnewred = 1

  sigma = 5. ; required significance level in contrast curves

  if (total(uint(filterflag)) gt 0.) then begin

    filtname = ['Wien', 'Catch', 'Gaussian Smoothing', 'Median Filtering', 'Mean Filtering', 'Gaussian Filtering', 'Fourier', 'Gaussian Fourier', 'Fourier Butterwoth']
    filtname_out = ['Wien', 'Ca', 'GaSm', 'MedF', 'MeF', 'GaF', 'Fo', 'GaFo', 'FoBu']

    ;for i=0,n_elements(filtname)-1 do print, strcompress(i+1), '  ', filtname[i]
    read, 'Select Filter: ', filtsel
;     filtsel = 4

;     if (n_elements(filterflag) eq 1) then filtmethod = [(filtname[filtsel-1])[0]]
;     if (n_elements(filterflag) eq 2 and filterflag[0] eq '0') then filtmethod = ['',(filtname[filtsel-1])[0]]
;     if (n_elements(filterflag) eq 2 and filterflag[0] eq '1') then filtmethod = [(filtname[filtsel-1])[0],'']

    if (n_elements(filterflag) eq 1) then begin
      filtmethod = [(filtname[filtsel-1])[0]]
      filtname_out = [(filtname_out[filtsel-1])[0]]
    endif
    if (n_elements(filterflag) eq 2 and filterflag[0] eq '0') then begin
      filtmethod = ['',(filtname[filtsel-1])[0]]
      filtname_out = ['',(filtname_out[filtsel-1])[0]]
    endif
    if (n_elements(filterflag) eq 2 and filterflag[0] eq '1') then begin
      filtmethod = [(filtname[filtsel-1])[0],'']
      filtname_out = [(filtname_out[filtsel-1])[0],'']
    endif

  endif else begin

    filtmethod = ['']

  endelse


  if (mth eq 0) then begin

    ;file selection

    if keyword_set(sphere) then begin

      if (flagnewred eq 0) then begin

	file = file_search(path+'*.fits', count=nfiles)

	out = strarr(nfiles)
	object = strarr(nfiles)
	type = strarr(nfiles)
	catg = strarr(nfiles)
	opti2 = strarr(nfiles)
	dit = strarr(nfiles)
	icor = strarr(nfiles)
	combind = strarr(nfiles)
	filt = strarr(nfiles)

	for i=0,nfiles-1 do begin

	  pos1 = strpos(file[i], '/', /reverse_search)
	  pos2 = strpos(file[i], '.fits', /reverse_search)
	  out[i] = strmid(file[i], pos1+1, pos2-pos1-1)

	  hdr = headfits(file[i], exten=0, /silent)

	  object[i] = strcompress(get_eso_keyword(hdr,'OBJECT'),/rem)
	  type[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DPR TYPE'),/rem)
	  catg[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DPR CATG'),/rem)
	  opti2[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS1 OPTI2 NAME'),/rem)
	  dit[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DET SEQ1 DIT'),/rem)
	  dit[i] = strmid(dit[i], 0, strpos(dit[i], '.'))
	  icor[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS COMB ICOR'),/rem)
	  combind[i] = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO INS4 FILT2 NAME'),/rem)	;Assembly for infrared neutral density
	  filt[i] = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO INS1 FILT NAME'),/rem)	;IRDIS filter unique ID

	endfor

	;selection of SCIENCE frames

	idx = where(catg eq 'SCIENCE')
	tfile = file[idx]
	tout = out[idx]
	ttype = type[idx]
	tcatg = catg[idx]
	tdit = dit[idx]
	topti2 = opti2[idx]
	ticor = icor[idx]
	tcombind = combind[idx]
	tfilt = filt[idx]

	print, ''
	for i=0,n_elements(idx)-1 do begin

	  print, strcompress(i+1,/rem)+' ', tout[i], tcatg[i], ttype[i], tdit[i], tfilt[i], topti2[i], tcombind[i], ticor[i], format='(a3, a32, a15, a17, a5, a8, a8, a10, a13)'
	  print, '----------------------------------------------------------------------------------------------------------------'

	endfor

	print, ''
	read, 'ONE Science Frame: ', sel
	scf = tfile[sel-1]
	read, 'ONE Object,Flux: ', sel
	fluxf = tfile[sel-1]

      endif else begin

	scf = datadir+'img_'+camera_filter[0]+'_dc.fits'
	fluxf = datadir+'PSF_'+camera_filter[0]+'.fits'

      endelse

      ;extract needed keywords

      hdrsc = headfits(scf, exten=0, /silent)
	sopti2 = strcompress(get_eso_keyword(hdrsc,'HIERARCH ESO INS1 OPTI2 NAME'),/rem)
	scombind = strcompress(get_eso_keyword(hdrsc, 'HIERARCH ESO INS4 FILT2 NAME'),/rem)	;Assembly for infrared neutral density
	sfilt = strcompress(get_eso_keyword(hdrsc, 'HIERARCH ESO INS1 FILT NAME'),/rem)	;IRDIS filter unique ID
	dit_sci = double(get_eso_keyword(hdrsc, 'HIERARCH ESO DET SEQ1 DIT'))
      hdrfl = headfits(fluxf, exten=0, /silent)
	fopti2 = strcompress(get_eso_keyword(hdrfl,'HIERARCH ESO INS1 OPTI2 NAME'),/rem)
	fcombind = strcompress(get_eso_keyword(hdrfl, 'HIERARCH ESO INS4 FILT2 NAME'),/rem)	;Assembly for infrared neutral density
	ffilt = strcompress(get_eso_keyword(hdrfl, 'HIERARCH ESO INS1 FILT NAME'),/rem)	;IRDIS filter unique ID
	dit_fl = double(get_eso_keyword(hdrfl, 'HIERARCH ESO DET SEQ1 DIT'))

      ;extract filter names for sph_ird_transmission.pro
      db_sci = sopti2
      bb_sci = sfilt
      if (bb_sci eq 'B_Ks') then bb_sci = 'B_K'
      if (scombind eq 'OPEN') then nd_sci = '0.0' $
	else nd_sci = strmid(scombind, strpos(scombind, '_', /reverse_search)+1, 3)

      db_fl = fopti2
      bb_fl = ffilt
      if (bb_fl eq 'B_Ks') then bb_fl = 'B_K'
      if (fcombind eq 'OPEN') then nd_fl = '0.0' $
	else nd_fl = strmid(fcombind, strpos(fcombind, '_', /reverse_search)+1, 3)

      transm_sci = 10.^sph_ird_transmission(bb_sci, db_sci, nd_sci)
      transm_fl = 10.^sph_ird_transmission(bb_fl, db_fl, nd_fl)

      ditfact = dit_sci/dit_fl

      if (bb_sci ne bb_fl and db_sci ne db_fl) then begin

	print, ''
	print, 'You screwed up your data.'
	stop

      endif

    endif

    ;========================================================================

    ;file selection

    if keyword_set(naco) then begin

      if (flagnewred eq 0) then begin

	file = file_search(path+'*.fits', count=nfiles)

	out = strarr(nfiles)
	object = strarr(nfiles)
	type = strarr(nfiles)
	catg = strarr(nfiles)
	opti1 = strarr(nfiles)
	opti3 = strarr(nfiles)
	opti6 = strarr(nfiles)
	opti7 = strarr(nfiles)
	dit = strarr(nfiles)
	naxis1 = strarr(nfiles)
	naxis2 = strarr(nfiles)

	for i=0,nfiles-1 do begin

	  pos1 = strpos(file[i], '/', /reverse_search)
	  pos2 = strpos(file[i], '.fits', /reverse_search)
	  out[i] = strmid(file[i], pos1+1, pos2-pos1-1)

	  hdr = headfits(file[i], exten=0, /silent)

	  object[i] = strcompress(get_eso_keyword(hdr,'OBJECT'),/rem)
	  naxis1[i] = strcompress(get_eso_keyword(hdr,'NAXIS1'),/rem)
	  naxis2[i] = strcompress(get_eso_keyword(hdr,'NAXIS2'),/rem)
	  type[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DPR TYPE'),/rem)
	  catg[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DPR CATG'),/rem)
	  dit[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DET DIT'),/rem)
	  dit[i] = sigfig(dit[i],2)
	  opti1[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS OPTI1 NAME'),/rem)	;mask
	  opti3[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS OPTI3 NAME'),/rem)	;pupil stop
	  opti6[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS OPTI6 NAME'),/rem)	;filter
	  opti7[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS OPTI7 NAME'),/rem)	;Objective

	endfor

	;selection of SCIENCE frames

	idx = where(catg eq 'SCIENCE')
	tfile = file[idx]
	tout = out[idx]
	ttype = type[idx]
	tcatg = catg[idx]
	tdit = dit[idx]
	topti1 = opti1[idx]
	topti3 = opti3[idx]
	topti6 = opti6[idx]
	topti7 = opti7[idx]
	tnaxis1 = naxis1[idx]
	tnaxis2 = naxis2[idx]

	print, ''
	for i=0,n_elements(idx)-1 do begin

	  print, strcompress(i+1,/rem), tout[i], tcatg[i], ttype[i], tdit[i], topti1[i], topti3[i], topti6[i], topti7[i], tnaxis1[i], tnaxis2[i], format='(a3, a40, a15, a17, a8, a8, a10, a10, 3a5)'
	  print, '------------------------------------------------------------------------------------------------------------------------------'

	endfor

	print, ''
	read, 'ONE Science Frame: ', sel
  ;     sel = 200
	scf = tfile[sel-1]
	read, 'ONE Flux Frame: ', sel
  ;     sel = 249
	fluxf = tfile[sel-1]

      endif else begin

	scf = datadir+'img_'+camera_filter[0]+'_dc.fits'
	fluxf = datadir+'PSF_'+camera_filter[0]+'.fits'

      endelse

      ;extract needed keywords

      hdrsc = headfits(scf, exten=0, /silent)
	dit_sci = double(get_eso_keyword(hdrsc, 'HIERARCH ESO DET DIT'))
	nd_sci = strcompress(get_eso_keyword(hdrsc, 'HIERARCH ESO INS OPTI3 ID'),/rem)
	mask = strcompress(get_eso_keyword(hdrsc,'HIERARCH ESO INS OPTI1 NAME'),/rem)
	if (mask eq 'AGPM') then agpm = 1 else agpm = 0
      hdrfl = headfits(fluxf, exten=0, /silent)
	dit_fl = double(get_eso_keyword(hdrfl, 'HIERARCH ESO DET DIT'))
	;dit_fl = 0.05
	nd_fl = strcompress(get_eso_keyword(hdrfl, 'HIERARCH ESO INS OPTI3 ID'),/rem)

      ;TO BE CLARIFIED: what is the influence of full_uszd?
;       if (nd_fl ne 'Full' and nd_fl ne 'ND_Long' and nd_fl ne 'Apo_165' and nd_fl ne 'Full_Uszd') then stop

      if (camera_filter[0] eq 'Lp') then begin

	if (nd_sci eq 'Full' or nd_sci eq 'Full_Uszd' and nd_fl eq 'ND_Long') then begin

	  filterpath = '/home/amueller/work/IDLlibs/AO/transmission_naco/Filters/'
	  readcol, filterpath+'L_prime.dat', wl, fl, format='d,d', /silent
	  readcol, filterpath+'ND_long.dat', wlnd, flnd, format='d,d', /silent

	  lambdainterp = dindgen(931.)/1000.+3.38d0
	  tl_interp = interpol(fl >0, wl, lambdainterp)
	  tnd_interp = interpol(flnd >0, wlnd, lambdainterp)

	  transm_fl = total(tl_interp)/total(tl_interp*tnd_interp)	;56.047

	  transm_sci = 1.
	  ;transm_fl = 50.	;see NACO manual p29

	endif else begin

	  transm_sci = 1.
	  transm_fl = 1.

	endelse
	
      endif
      
      if (camera_filter[0] eq 'Mp') then begin

	if (nd_sci eq 'Full' or nd_sci eq 'Full_Uszd' and nd_fl eq 'ND_Long') then begin

	  filterpath = '/home/amueller/work/IDLlibs/AO/transmission_naco/Filters/'
	  readcol, filterpath+'M_prime.dat', wl, fl, format='d,d', /silent
	  readcol, filterpath+'ND_long.dat', wlnd, flnd, format='d,d', /silent

	  lambdainterp = dindgen(931.)/1000.+4.35d0
	  tl_interp = interpol(fl >0, wl, lambdainterp)
	  tnd_interp = interpol(flnd >0, wlnd, lambdainterp)

	  transm_fl = total(tl_interp)/total(tl_interp*tnd_interp)

	  transm_sci = 1.

	endif else begin

	  transm_sci = 1.
	  transm_fl = 1.

	endelse
	
      endif

      if (camera_filter[0] eq 'Ks') then begin

	if (nd_sci eq 'Full' and nd_fl eq 'ND_Short') then begin

	  filterpath = '/home/amueller/work/IDLlibs/AO/transmission_naco/Filters/'
	  readcol, filterpath+'Ks_CONICA.dat', wl, fl, format='d,d', /silent
	  readcol, filterpath+'ND_short.dat', wlnd, flnd, format='d,d', /silent

	  lambdainterp = dindgen(561.)/1000.+1.88d0
	  tl_interp = interpol(fl >0, wl, lambdainterp)
	  tnd_interp = interpol(flnd >0, wlnd, lambdainterp)

	  transm_fl = total(tl_interp)/total(tl_interp*tnd_interp)

	  transm_sci = 1.

	endif else begin

	  transm_sci = 1.
	  transm_fl = 1.

	endelse
	
      endif

      ditfact = dit_sci/dit_fl

    endif

    ;========================================================================

    ;file selection

    if keyword_set(nici) then begin

      scf = datadir+'img_'+camera_filter[0]+'_dc.fits'
      fluxf = datadir+'PSF_'+camera_filter[0]+'.fits'

      hdrsc = headfits(scf, exten=0, /silent)
      	dit_sci = double(get_eso_keyword(hdrsc, 'ITIME_R'))
	dit_fl = dit_sci

      ;DUMMY
      ditfact = 1.
      ;dit_sci = 1.
      ;dit_fl = 1.
      transm_sci = 1.
      transm_fl = 1.

    endif
    
    ;========================================================================

    ;file selection

    if keyword_set(nirc2) then begin

      scf = datadir+'img_'+camera_filter[0]+'_dc.fits'
      fluxf = datadir+'PSF_'+camera_filter[0]+'.fits'

      hdrsc = headfits(scf, exten=0, /silent)
      	dit_sci = double(get_eso_keyword(hdrsc, 'ITIME'))
	dit_fl = dit_sci

      ;DUMMY
      ditfact = 1.
      ;dit_sci = 1.
      ;dit_fl = 1.
      transm_sci = 1.
      transm_fl = 1.

    endif
    
    ;========================================================================

    ;file selection

    if keyword_set(magao) then begin

      scf = datadir+'img_'+camera_filter[0]+'_dc.fits'
      fluxf = datadir+'PSF_'+camera_filter[0]+'.fits'

      hdrsc = headfits(scf, exten=0, /silent)
      	dit_sci = double(get_eso_keyword(hdrsc, 'DIT'))
	dit_fl = dit_sci

      ;DUMMY
      ditfact = 1.
      ;dit_sci = 1.
      ;dit_fl = 1.
      transm_sci = [1.,1.]
      transm_fl = [1.,1.]

    endif

    ;========================================================================

    ;file selection

    if keyword_set(lmir) then begin

      scf = datadir+'img_'+camera_filter[0]+'_dc.fits'
      fluxf = datadir+'PSF_'+camera_filter[0]+'.fits'

      hdrsc = headfits(scf, exten=0, /silent)
      	dit_sci = double(get_eso_keyword(hdrsc, 'ITIME'))
      hdrfl = headfits(fluxf, exten=0, /silent)
	dit_fl = double(get_eso_keyword(hdrfl, 'ITIME'))
	ffilt = strcompress(get_eso_keyword(hdrfl, 'LMIR_FW2'))

      ;assuming that science was taken w/o ND filter
      ;assuming that ND filter has wavelength independent troughput of 10%
      if (ffilt eq 'ND1.0-T10') then transm_fl = [9.03,9.03] else transm_fl = [1.,1.]

      ditfact = dit_sci/dit_fl
      transm_sci = [1.,1.]

    endif

  endif

  ;========================================================================

  if not keyword_set(plotonly) then begin

    if (mth eq 0) then begin

      qcut = ''
      if (fcpflag eq '1') then begin

	if keyword_set(cut) then qcut = 'y' else qcut = 'n'
	;read, 'Cut out CC for computation of contrast curve (y/n): ', qcut

      endif else begin

	qcut = 'n'

      endelse

    endif

  endif

  ;========================================================================

  ;obsolete
  ; normal = 8000 ; normalization by photometric variations (TBC)
  ; normal = 1
  ; box = 1

  ;========================================================================

  procdir = datadir

  tg_name = camera_filter
  if (keyword_set(sphere) and n_elements(camera_filter) eq 2) then filepsf = ['PSF_'+camera_filter[0]+'.fits', 'PSF_'+camera_filter[1]+'.fits'] else filepsf = ['PSF_'+camera_filter[0]+'.fits']
  if keyword_set(naco) then filepsf = ['PSF_'+camera_filter[0]+'.fits']
  if keyword_set(nici) then filepsf = ['PSF_'+camera_filter[0]+'.fits']
	if keyword_set(nirc2) then filepsf = ['PSF_'+camera_filter[0]+'.fits']
  if keyword_set(lmir) then filepsf = ['PSF_'+camera_filter[0]+'.fits', 'PSF_'+camera_filter[1]+'.fits']
  if keyword_set(magao) then filepsf = ['PSF_'+camera_filter[0]+'.fits', 'PSF_'+camera_filter[1]+'.fits']
  
  ;========================================================================

  ;read in to get dimensions
  tmp = readfits(datadir+'img_'+tg_name[0]+'_dc.fits', 0, /silent)
  sz = size(tmp)

  dim_init = sz[1] ; Input image size (not used, unless WCS data are considered)
  cx = dim_init/2.
  cy = dim_init/2.

  ;some simple checks if parameters are ok
  ; if (method eq 'c') then begin
  ; 
  ;   if (max(zone) gt dim) then stop	;zone cannot be larger than picture
  ;   if (max(klip) gt sz[3]) then stop	;klip cannot be larger than number of frames
  ; 
  ; endif

  ;========================================================================

  if keyword_set(sphere) then plsc = 0.01225d0 ; plate scale [arcsec/px]
  if (keyword_set(naco) and camera_filter[0] eq 'Lp' or camera_filter[0] eq 'Mp') then plsc = 0.02719d0	;L27
  if (keyword_set(naco) and camera_filter[0] eq 'Ks') then plsc = 0.013221d0	;S13
  diam = 8.2 ; telescope diameter [m]
  if keyword_set(nici) then begin
    plsc = 0.017958d0	;L
    diam = 8.1
  endif
  if keyword_set(lmir) then begin
    plsc = 0.010707d0	;L
    diam = 8.4
  endif
  if keyword_set(magao) then begin
    plsc = 0.007851d0	;Vis
    diam = 6.5
  endif
  if keyword_set(nirc2) then begin
		plsc = 0.009952d0
		diam = 10.
  endif


  ;========================================================================

  ndelta = n_elements(delta)
  if (method eq 'a' or method eq 'c') then ndelta = 1

  nklip = n_elements(klip)

  if (method eq 'c') then nzone = n_elements(zone) else nzone = n_elements(step_init)
  if (method eq 'p' or method eq 'v' or method eq 'vc' or method eq 'loci') then nzone = n_elements(step_init)

  if (method eq 'a' or method eq 'sa') then begin

    nklip = 1
    nzone = 1

  endif

  nrin = n_elements(rin_init)

  if (method eq 'p' or method eq 'v' or method eq 'vc' or method eq 'c') then rootresdir = datadir+'PCA_Reduction/'
  if (method eq 'l') then rootresdir = datadir+'LLSG_Reduction/'
  if (method eq 'loci') then rootresdir = datadir+'LOCI_Reduction/'
  if (method eq 'a') then rootresdir = datadir+'ADI_Reduction/'
  if (method eq 'sa') then rootresdir = datadir+'sADI_Reduction/'

  ;========================================================================

  spawn, 'mkdir -p '+rootresdir

  for j=0,nklip-1 do begin

    for k=0,nzone-1 do begin

      for l=0,nrin-1 do begin

	for m=0,ndelta-1 do begin

	  for i=0,n_elements(filterflag)-1 do begin

	    for xx=0,n_elements(camera_filter)-1 do begin


	      tg_name_ori = tg_name[xx]

	      ;H23 not the same as H32
	      if (camera_filter[xx] eq 'K1' and xx eq 0) then lambda = 2.1025d-6
	      if (camera_filter[xx] eq 'H2' and xx eq 0) then lambda = 1.5888d-6
	      if (camera_filter[xx] eq 'H3' and xx eq 0) then lambda = 1.6653d-6
	      if (camera_filter[xx] eq 'J2' and xx eq 0) then lambda = 1.1895d-6

	      if (camera_filter[xx] eq 'K2' and xx eq 1) then lambda = 2.2550d-6
	      if (camera_filter[xx] eq 'H2' and xx eq 1) then lambda = 1.5890d-6
	      if (camera_filter[xx] eq 'H3' and xx eq 1) then lambda = 1.6671d-6
	      if (camera_filter[xx] eq 'J3' and xx eq 1) then lambda = 1.2698d-6

	      if (camera_filter[xx] eq 'YL' and xx eq 0) then lambda = 1.0425d-6
	      if (camera_filter[xx] eq 'YR' and xx eq 0) then lambda = 1.0425d-6
	      if (camera_filter[xx] eq 'JL' and xx eq 0) then lambda = 1.2575d-6
	      if (camera_filter[xx] eq 'JR' and xx eq 0) then lambda = 1.2575d-6
	      if (camera_filter[xx] eq 'HL' and xx eq 0) then lambda = 1.6255d-6
	      if (camera_filter[xx] eq 'HR' and xx eq 0) then lambda = 1.6255d-6
	      if (camera_filter[xx] eq 'KL' and xx eq 0) then lambda = 2.1813d-6
	      if (camera_filter[xx] eq 'KR' and xx eq 0) then lambda = 2.1813d-6
	      
	      if (camera_filter[xx] eq 'PaB') then lambda = 1.283d-6

	      if (camera_filter[xx] eq 'Lp') then lambda = 3.80d-6
	      if (camera_filter[xx] eq 'Mp') then lambda = 4.78d-6
	      if (camera_filter[xx] eq 'Ks') then lambda = 2.20d-6

	      if (camera_filter[xx] eq 'dx' or camera_filter[xx] eq 'sx') then lambda = 3.70d-6
	      
	      if (camera_filter[xx] eq 'Ha') then lambda = 656.d-9
	      if (camera_filter[xx] eq 'CntHa') then lambda = 643.d-9
	      
	      if (keyword_set(nirc2) and camera_filter[xx] eq 'Lp') then lambda = 3.776d-6
	      if (keyword_set(nirc2) and camera_filter[xx] eq 'Ms') then lambda = 4.67050d-6
	      
	      fwhm = lambda/diam*206265.d0/plsc ; fhwm of the psf in px, assuming no Lyot stop
						;The number 206,265 is approximately equal to the number of arcseconds in a circle (1,296,000), divided by 2π.  360.*60.*60/(2.*!DPI)
	      highpass = 2.d0*fwhm ; parameter of high pass filter
	      lowpass = fwhm/2.	;if parameters change, manually change it in fcp_am

	      ;========================================================================

	      procdir = datadir
	      tg_name_dyn = tg_name_ori
	      tg_name_bas = tg_name_ori

	      if not keyword_set(plotonly) then begin

		;========================================================================

		MERGING_WINDOW_V5, cx, cy, cx, cy, dim, /resize_only

		;========================================================================

		;read in PSF and scale it

		psf_orig = mrdfits(datadir+'PSF_'+tg_name_bas+'.fits', 0, hdrpsf, /silent)
		psf = psf_orig*ditfact*transm_fl[xx]/transm_sci[xx]
		writefits, procdir+'img_'+tg_name_bas+'_scaledPSF.fits', psf, hdrpsf

		circint_MOD, psf, (n_elements(psf[*,0])-1.)/2, (n_elements(psf[0,*])-1.)/2, fwhm, tot, mtot, meantot, maxpx, sdtot, npx, totdif, npxdif, t8;
		psfflux = tot[0]
		meanpsfflux = meantot[0]

                if keyword_set(naco) then begin

                  filepsfall = file_search(path+'../Reduced/'+'PSF_2*.fits', count=npsf)
                  ;psfall = dblarr(n_elements(psf[*,0]), n_elements(psf[*,0]), npsf)
		  tmphdr = headfits(filepsfall[0], exten=0)
		  tmpdim = get_eso_keyword(tmphdr, 'NAXIS1')
		  psfall = dblarr(tmpdim, tmpdim, npsf)

                  for pp=0,npsf-1 do begin

                   psfall[*,*,pp] = mrdfits(filepsfall[pp],0,/silent)
                   psfall[*,*,pp] = psfall[*,*,pp]*ditfact*transm_fl[xx]/transm_sci[xx]

                  endfor

                  writefits, procdir+'img_'+tg_name_bas+'_scaledPSF_frames.fits', psfall, hdrpsf

                endif

                if keyword_set(nici) then begin

                  filepsfall = file_search(path+'../Reduced/'+'PSF_2*.fits', count=npsf)
                  ;psfall = dblarr(n_elements(psf[*,0]), n_elements(psf[*,0]), npsf)
		  tmphdr = headfits(filepsfall[0], exten=0)
		  tmpdim = get_eso_keyword(tmphdr, 'NAXIS1')
		  psfall = dblarr(tmpdim, tmpdim, npsf)

                  for pp=0,npsf-1 do begin

                   psfall[*,*,pp] = mrdfits(filepsfall[pp],0,/silent)
                   psfall[*,*,pp] = psfall[*,*,pp]*ditfact*transm_fl[xx]/transm_sci[xx]

                  endfor

                  writefits, procdir+'img_'+tg_name_bas+'_scaledPSF_frames.fits', psfall, hdrpsf

                endif
                
                if keyword_set(nirc2) then begin

                  filepsfall = file_search(path+'../Reduced/'+'PSF_2*.fits', count=npsf)
                  ;psfall = dblarr(n_elements(psf[*,0]), n_elements(psf[*,0]), npsf)
		  tmphdr = headfits(filepsfall[0], exten=0)
		  tmpdim = get_eso_keyword(tmphdr, 'NAXIS1')
		  psfall = dblarr(tmpdim, tmpdim, npsf)

                  for pp=0,npsf-1 do begin

                   psfall[*,*,pp] = mrdfits(filepsfall[pp],0,/silent)
                   psfall[*,*,pp] = psfall[*,*,pp]*ditfact*transm_fl[xx]/transm_sci[xx]

                  endfor

                  writefits, procdir+'img_'+tg_name_bas+'_scaledPSF_frames.fits', psfall, hdrpsf

                endif

                if keyword_set(lmir) then begin

                  filepsfall = file_search(path+'../Reduced/'+'PSF_'+camera_filter[xx]+'_*.fits', count=npsf)
		  tmphdr = headfits(filepsfall[0], exten=0)
		  tmpdim = get_eso_keyword(tmphdr, 'NAXIS1')
		  psfall = dblarr(tmpdim, tmpdim, npsf)

                  for pp=0,npsf-1 do begin

                   psfall[*,*,pp] = mrdfits(filepsfall[pp],0,/silent)
                   psfall[*,*,pp] = psfall[*,*,pp]*ditfact*transm_fl[xx]/transm_sci[xx]

                  endfor

                  writefits, procdir+'img_'+tg_name_bas+'_scaledPSF_frames.fits', psfall, hdrpsf

                endif
                
                if keyword_set(sphere) then begin

                  filepsfall = file_search(path+'PSF_'+camera_filter[xx]+'_*.fits', count=npsf)

		  naxis3 = dblarr(npsf)
		  for pp=0,npsf-1 do begin

		    hdr = headfits(filepsfall[pp], exten=0, /silent)
		    naxis3[pp] = get_eso_keyword(hdr, 'NAXIS3')

		  endfor

		  if (naxis3[0] eq 0.) then naxis3 = 1

		  ;if (naxis3 eq 0.) then naxis3 = 1
		  ;psfall = dblarr(n_elements(psf[*,0]), n_elements(psf[*,0]), total(naxis3))
		  tmphdr = headfits(filepsfall[0], exten=0)
		  tmpdim = get_eso_keyword(tmphdr, 'NAXIS1')
		  psfall = dblarr(tmpdim, tmpdim, total(naxis3))

		  for pp=0,npsf-1 do begin

		    if (pp eq 0) then psfall[*,*,0:naxis3[pp]-1] = mrdfits(filepsfall[pp], 0, /silent)
		    if (pp gt 0) then psfall[*,*,total(naxis3[0:pp-1]):total(naxis3[0:pp-1])+naxis3[pp]-1] = mrdfits(filepsfall[pp], 0, /silent)

		  endfor

		  for pp=0,n_elements(psfall[0,0,*])-1 do psfall[*,*,pp] = psfall[*,*,pp]*ditfact*transm_fl[xx]/transm_sci[xx]

                  writefits, procdir+'img_'+tg_name_bas+'_scaledPSF_frames.fits', psfall, hdrpsf

                endif


		;scale image cube
		objtmp = mrdfits(procdir+'img_'+tg_name_dyn+'_dc.fits', 0, header, /silent)
		objtmp = objtmp/max(psf);   meanpsfflux; max(psf)	;psfflux

		if (selins eq 'n') then begin

		  if (agpm eq 1) then begin

		    r_fwhm = dindgen(dim)/fwhm
		    ;equation from D. Mawet, reference Delacroix2013, A&A 553, A98
		    ;r_fwhm = r_fwhm*0.85 ; correct the effective fwhm
		    ;coro_thr = 1.d0-exp(-0.831*(r_fwhm)^2) ; transmission profile

		    ;measured by ISPY
		    coro_thr = 1.d0-exp(-0.58566*(r_fwhm)^1.32979) ; transmission profile

		    coro_thr_2d = dblarr(dim,dim)

		    for xct=0d,dim-1 do begin

		      for yct=0d,dim-1 do begin

			rect_coord = [xct-dim/2, yct-dim/2]
			polar_coord = cv_coord(FROM_RECT=rect_coord, /TO_POLAR, /DEGREES)
			if abs(polar_coord[1]) gt dim/2-1 then coro_thr_2d[xct,yct] = 0 else coro_thr_2d[xct,yct] = coro_thr[abs(polar_coord[1])]

		      endfor

		    endfor

		    coro_thr_2d[where(coro_thr_2d eq 0.)] = 1.
		    for xy=0,n_elements(objtmp[0,0,*])-1 do objtmp[*,*,xy] = objtmp[*,*,xy]/coro_thr_2d

		  endif

		endif


		tg_name_dyn = tg_name_dyn+'_scaled'
		writefits, procdir+'img_'+tg_name_dyn+'_dc.fits', objtmp, header

		;========================================================================

		if (recenter eq 'y') then begin

		  im = readfits(procdir+'img_'+tg_name_dyn+'_dc.fits', header, /silent)

		  ;define radii, keep it SINGLE radii for the moment
		  if (rin_init[l] ne 0) then begin

		    mask_t = shift(dist(dim),dim/2,dim/2)
		    mask = mask_t ge fwhm

		    ;rin = [ceil(fwhm*rin_init[l]),22,60]
		    ;rout = [22,60,dim/2-5]
		    ;rin = [ceil(fwhm*rin_init[l])]
		    ;rin = [20.]
		    ;rout = [dim/2-5]	;[90]
		    if keyword_set(sphere) then begin

		      rin = [50.]
		      rout = [90.]

		    endif

		    if keyword_set(naco) then begin

		      if (agpm eq 1) then begin

			rin = [10.]
			rout = [30.]

		      endif else begin

			rin = [8.]
			rout = [18.]

		      endelse

		    endif


		  endif else begin

		    ;rin = [10,22,60]
		    ;rout = [22,60,dim/2-5]
		    rin = [20]
		    rout = [dim/2-5]	;[90]

		    if keyword_set(sphere) then begin

		      rin = [50.]
		      rout = [100.]

		    endif

		    if keyword_set(naco) then begin

		      rin = [10.]
		      rout = [30.]

		    endif

		  endelse

		  psfref = median(im, dim=3, /even)
		  diff = dblarr(n_elements(rin),dim,dim,(size(im))[3])

		  ;only translation and intensity, no rotation or magnification
		  subxy = dblarr(2,(size(im))[3])
		  for ii=0,(size(im))[3]-1 do begin

		    ;mode=0: pr, sr
		    ;mode=1: pr, sr, ar
		    ;mode=2: pr, sr, ar, zr
		    ; pr : is the translation range (roughly) in pixels.
		    ; sr : is the instensity scaling
		    ; ar : is the rotation range in degrees.
		    ; zr : is the zooming range
		    ; ftol : is the toleration to stop the fitting.

		    mode = 0
		    w = get_w(psfref, rin, rout)
		    for jj=0,n_elements(rin)-1 do diff[jj,*,*,ii] = transmatch9(psfref, im[*,*,ii], w, dim/2., dim/2., mode=mode, pr=2.,sr=0.1,ar=2.,zr=0.05, rin=rin[jj], rout=rout[jj], ftol=1.e-6,/nodisp,/silent)

		    rcflag = '_Recenter'+strcompress(mode,/rem)

		  endfor

		  diff = reform(diff)	;only using 1 radii for the moment

		  if (method ne 'a' or method ne 'sa') then begin

		    for ii=0,(size(im))[3]-1 do im[*,*,ii] = psfref-diff[*,*,ii]

		  endif else begin

		    for ii=0,(size(im))[3]-1 do im[*,*,ii] = psfref-diff[*,*,ii]
		    diff2 = dblarr(dim, dim, (size(im))[3])
		    for ii=0,(size(im))[3]-1 do diff2[*,*,ii] = im[*,*,ii]-psfref
		    im = diff2

		  endelse

		  tg_name_dyn = tg_name_dyn+'_shift'
		  writefits, datadir+'img_'+tg_name_dyn+'_dc.fits', im, header

		endif else begin

		  rcflag = ''

		endelse

		;========================================================================

    ; 	    if (filterflag[i] eq '1') then MFILTER_V4b, HIGHPASS=cutoff_l, LOWPASS=fwhm/2, DISPLAY=display, /CA;/gauss;, /CA;, /MEDIAN_F
		if (filterflag[i] eq '1') then MFILTER_V4b, filtmethod[i], HIGHPASS=highpass, LOWPASS=lowpass, DISPLAY=display;, /CA;/gauss;, /CA;, /MEDIAN_F

; 		if (filterflag[i] eq '1') then MFILTER_V4b, filtmethod[i], HIGHPASS=17, DISPLAY=display;, /CA;/gauss;, /CA;, /MEDIAN_F

	      endif

	      ;========================================================================

	      ;define output directories

	      if (filterflag[i] eq '0' and method eq 'c') then outdir = 'PCA_noFilter_KLIP'+strcompress(klip[j],/rem)+'_Zone'+strcompress(zone[k],/rem)+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)
	      if (filterflag[i] eq '1' and method eq 'c') then outdir = 'PCA_Filter_KLIP'+strcompress(klip[j],/rem)+'_Zone'+strcompress(zone[k],/rem)+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)

	      if (pointdisk eq 'd') then begin

		if (filterflag[i] eq '0' and method eq 'v') then outdir = 'sPCAvar_noFilter_KLIP'+strcompress(klip[j],/rem)+'_Step'+strcompress(uint(step_init[k]),/rem)+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+'_Delta'+strcompress(sigfig(delta[m],2),/rem)+rcflag
		if (filterflag[i] eq '1' and method eq 'v') then outdir = 'sPCAvar_'+filtname_out[i]+'Filter_KLIP'+strcompress(klip[j],/rem)+'_Step'+strcompress(uint(step_init[k]),/rem)+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+'_Delta'+strcompress(sigfig(delta[m],2),/rem)+rcflag

		if (filterflag[i] eq '0' and method eq 'vc') then outdir = 'sPCAvarcorr_noFilter_KLIP'+strcompress(klip[j],/rem)+'_Step'+strcompress(uint(step_init[k]),/rem)+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+'_Delta'+strcompress(sigfig(delta[m],2),/rem)+rcflag
		if (filterflag[i] eq '1' and method eq 'vc') then outdir = 'sPCAvarcorr_'+filtname_out[i]+'Filter_KLIP'+strcompress(klip[j],/rem)+'_Step'+strcompress(uint(step_init[k]),/rem)+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+'_Delta'+strcompress(sigfig(delta[m],2),/rem)+rcflag

	      endif else begin

		if (filterflag[i] eq '0' and method eq 'v') then outdir = 'sPCAvar_noFilter_KLIP'+strcompress(klip[j],/rem)+'_Step'+strcompress(uint(step_init[k]),/rem)+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+'_Delta'+strcompress(sigfig(delta[m],2),/rem)+rcflag
		if (filterflag[i] eq '1' and method eq 'v') then outdir = 'sPCAvar_'+filtname_out[i]+'Filter_KLIP'+strcompress(klip[j],/rem)+'_Step'+strcompress(uint(step_init[k]),/rem)+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+'_Delta'+strcompress(sigfig(delta[m],2),/rem)+rcflag

		if (filterflag[i] eq '0' and method eq 'vc') then outdir = 'sPCAvarcorr_noFilter_KLIP'+strcompress(klip[j],/rem)+'_Step'+strcompress(uint(step_init[k]),/rem)+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+'_Delta'+strcompress(sigfig(delta[m],2),/rem)+rcflag
		if (filterflag[i] eq '1' and method eq 'vc') then outdir = 'sPCAvarcorr_'+filtname_out[i]+'Filter_KLIP'+strcompress(klip[j],/rem)+'_Step'+strcompress(uint(step_init[k]),/rem)+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+'_Delta'+strcompress(sigfig(delta[m],2),/rem)+rcflag

	      endelse

	      if (filterflag[i] eq '0' and method eq 'p' and smethod eq '3') then outdir = 'sPCAv13_noFilter_KLIP'+strcompress(klip[j],/rem)+'_Step'+strcompress(uint(step_init[k]),/rem)+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+'_Delta'+strcompress(sigfig(delta[m],2),/rem)+rcflag
	      if (filterflag[i] eq '1' and method eq 'p' and smethod eq '3') then outdir = 'sPCAv13_'+filtname_out[i]+'Filter_KLIP'+strcompress(klip[j],/rem)+'_Step'+strcompress(uint(step_init[k]),/rem)+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+'_Delta'+strcompress(sigfig(delta[m],2),/rem)+rcflag

	      if (filterflag[i] eq '0' and method eq 'l') then outdir = 'LLSG_noFilter_Rank'+strcompress(uint(klip[j]),/rem)+'_Thresh'+strcompress(uint(thresh),/rem)+rcflag
	      if (filterflag[i] eq '1' and method eq 'l') then outdir = 'LLSG_'+filtname_out[i]+'Filter_Rank'+strcompress(uint(klip[j]),/rem)+'_Thresh'+strcompress(uint(thresh),/rem)+rcflag

	      if (filterflag[i] eq '0' and method eq 'p' and smethod eq '5') then outdir = 'sPCAv15_noFilter_KLIP'+strcompress(klip[j],/rem)+'_Step'+strcompress(uint(step_init[k]),/rem)+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+'_Delta'+strcompress(sigfig(delta[m],2),/rem)+rcflag
	      if (filterflag[i] eq '1' and method eq 'p' and smethod eq '5') then outdir = 'sPCAv15_'+filtname_out[i]+'Filter_KLIP'+strcompress(klip[j],/rem)+'_Step'+strcompress(uint(step_init[k]),/rem)+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+'_Delta'+strcompress(sigfig(delta[m],2),/rem)+rcflag

	      if (filterflag[i] eq '0' and method eq 'loci') then outdir = 'LOCI_noFilter_KLIP'+strcompress(klip[j],/rem)+'_Step'+strcompress(uint(step_init[k]),/rem)+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+'_Delta'+strcompress(sigfig(delta[m],2),/rem)+'_Na'+strcompress(uint(Na),/rem)+rcflag
	      if (filterflag[i] eq '1' and method eq 'loci') then outdir = 'LOCI_'+filtname_out[i]+'Filter_KLIP'+strcompress(klip[j],/rem)+'_Step'+strcompress(uint(step_init[k]),/rem)+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+'_Delta'+strcompress(sigfig(delta[m],2),/rem)+'_Na'+strcompress(uint(Na),/rem)+rcflag

	      if (filterflag[i] eq '0' and method eq 'a') then outdir = 'ADI_noFilter'+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+rcflag
	      if (filterflag[i] eq '1' and method eq 'a') then outdir = 'ADI_'+filtname_out[i]+'Filter'+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+rcflag

	      if (filterflag[i] eq '0' and method eq 'sa') then outdir = 'sADI_noFilter'+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+'_Delta'+strcompress(sigfig(delta[m],2),/rem)+'_D'+sigfig(distdelta,2)+rcflag
	      if (filterflag[i] eq '1' and method eq 'sa') then outdir = 'sADI_'+filtname_out[i]+'Filter'+'_Rin'+strcompress(sigfig(rin_init[l],2),/rem)+'_Delta'+strcompress(sigfig(delta[m],2),/rem)+'_D'+sigfig(distdelta,2)+rcflag



	      resdir = rootresdir+outdir+'/'
	      spawn, 'mkdir -p '+resdir

	      if not keyword_set(plotonly) then begin

; 		if (method eq 'p') then begin
; 
; 		  resdirev1 = resdir+'Eigenvalues/'
; 		  spawn, 'mkdir -p '+resdirev1
; 		  resdirev = resdirev1+camera_filter[xx]+'/'
; 		  spawn, 'mkdir -p '+resdirev
; 
; 		endif
		if (method eq 'v' or method eq 'vc') then begin

		  resdirkd1 = resdir+'KLIP_Delta/'
		  spawn, 'mkdir -p '+resdirkd1
		  resdirkd = resdirkd1+camera_filter[xx]+'/'
		  spawn, 'mkdir -p '+resdirkd

		endif

		;========================================================================

		;Reduction

		if (method eq 'c') then reduced_img_median = PCA_ADI_V39_AM(cx, cy, dim, dim_init, fwhm, klip[j], zone[k], rin_init[l], delta[m], right_handed, 0, 0, DISPLAY=display, /verbose, NORMAL=normal)
		if (method eq 'v') then reduced_img_median = SPCA_ADI_VAR(cx, cy, dim, dim_init, fwhm, plsc, klip[j], rin_init[l], step_init[k], delta[m], right_handed, 0, 0, pointdisk, 0, 0, /verbose, /output, section=section)
		if (method eq 'vc') then reduced_img_median = SPCA_ADI_VAR_CORR(cx, cy, dim, dim_init, fwhm, plsc, klip[j], rin_init[l], step_init[k], delta[m], right_handed, 0, 0, pointdisk, 0, 0, /verbose, /output, section=section)
		if (method eq 'p' and smethod eq '3') then reduced_img_median = SPCA_ADI_V13(cx, cy, dim, dim_init, fwhm, plsc, klip[j], rin_init[l], step_init[k], delta[m], right_handed, 0, 0, NORMAL=normal, DISPLAY=display, /verbose)
		if (method eq 'p' and smethod eq '5') then reduced_img_median = SPCA_ADI_V15(cx, cy, dim, dim_init, fwhm, plsc, klip[j], rin_init[l], step_init[k], delta[m], right_handed, 0, 0, NORMAL=normal, DISPLAY=display, /verbose)
		if (method eq 'l') then reduced_img_median = llsg_am(cx, cy, dim, dim_init, fwhm, plsc, right_handed, klip[j], thresh, max_iter, low_rank_mode, thresh_mode, rin_init[l], delta[m], 0, 0)
		if (method eq 'loci') then reduced_img_median = loci_adi_v19(cx, cy, dim, dim_init, plsc, fwhm, rin_init[l], step_init[k], klip[j], Na, g, delta[m], right_handed, 0, 0, NORMAL=normal, /verbose, /svd, /mask_zone_s);, DISPLAY=display
		;/svd, /nnls, /bvls, /damped)
		if (method eq 'a') then reduced_img_median = adi_v01(cx, cy, dim, dim_init, plsc, fwhm, rin_init[l], right_handed, recenter, 0, 0, /verbose)
		if (method eq 'sa') then reduced_img_median = sadi(cx, cy, dim, dim_init, plsc, fwhm, rin_init[l], delta[m], distdelta, right_handed, recenter, 0, 0, /verbose)

		;========================================================================

		;convolve reduced image with PSF
		;psf = mrdfits(filepsf, 0)
		; cimg = convol(reduced_img_median, psf, /center, /edge_truncate)
		cimg = filter_image(reduced_img_median, fwhm_gaussian=0.5*fwhm);0.5*	;like pynpoint
		fn = 'img_'+tg_name_dyn+'_median_convolved.fits'
		header = headfits(procdir+'img_'+tg_name_ori+'_dc.fits', exten=0, /silent)
		mask_t = shift(dist(dim), dim/2., dim/2.)
		mask = mask_t ge fwhm*rin_init[l]
		writefits, datadir+fn, cimg*mask, header

		;========================================================================

		if (canddet eq '1' and fcpflag eq '1') then candidate_detection, star, camera_filter[xx], filterflag[i], dim, plsc, fwhm

		;========================================================================

		print, ''
		print, ''

		if (fcpflag eq '1') then begin

		  ;noise curve (no normalization to alpha_r)
		  ;radin = (rin_init[l]+0.5 > 2.0)*fwhm	;px

		  if keyword_set(sphere) then noise_curve_am, reduced_img_median, dim/2, dim/2, fwhm, plsc, dim, sigma, rin_init[l], qcut, DISPLAY=display;, box=1;, /GAUSSFILT
		  if keyword_set(naco) then noise_curve_am, reduced_img_median, dim/2, dim/2, fwhm, plsc, dim, sigma, rin_init[l], qcut, DISPLAY=display;, box=1;, /GAUSSFILT, agpm=agpm
		  if keyword_set(nici) then noise_curve_am, reduced_img_median, dim/2, dim/2, fwhm, plsc, dim, sigma, rin_init[l], qcut, DISPLAY=display;, box=1;, /GAUSSFILT, agpm=agpm
		  if keyword_set(lmir) then noise_curve_am, reduced_img_median, dim/2, dim/2, fwhm, plsc, dim, sigma, rin_init[l], qcut, DISPLAY=display;, box=1;, /GAUSSFILT, agpm=agpm

		endif

		;========================================================================

		;write out log file

		openw, lun, datadir+'setup_'+camera_filter[xx]+'.txt', width=2000, /get_lun

		  printf, lun, 'Target: '+star
		  printf, lun, ''
		  printf, lun, 'Filter: '+camera_filter[xx]
		  printf, lun, 'Wavelength/um: '+sigfig(lambda*1.d6,5)
		  printf, lun, 'PlateScale/mas: '+sigfig(plsc*1.d3,4)
		  printf, lun, 'FWHM_PSF/px: '+strcompress(fwhm,/rem)
		  printf, lun, ''
		  printf, lun, 'Star_xpos/px(IDL): '+sigfig(cx,4)
		  printf, lun, 'Star_ypos/px(IDL): '+sigfig(cy,4)
		  printf, lun, 'Dimension/px: '+strcompress(dim,/rem)
		  printf, lun, ''
		  if (method eq 'c') then printf, lun, 'Algorithm: '+'PCA orig'
		  if (method eq 'v') then printf, lun, 'Algorithm: '+'sPCA Var'
		  if (method eq 'vc') then printf, lun, 'Algorithm: '+'sPCA Var Corr'
		  if (method eq 'p') then printf, lun, 'Algorithm: '+'sPCA'
		  if (method eq 'l') then printf, lun, 'Algorithm: '+'LLSG'
		  if (method eq 'loci') then printf, lun, 'Algorithm: '+'LOCI' 
		  if (method eq 'a') then printf, lun, 'Algorithm: '+'ADI'
		  if (method eq 'sa') then printf, lun, 'Algorithm: '+'sADI'
		  printf, lun, ''
		  if (recenter eq 'y') then printf, lun, 'Recentering: yes' else printf, lun, 'Recentering: no'
		  if (filterflag[i] eq '0') then begin
		    printf, lun, 'Spatial Filter: '+'None'
		    ;printf, lun, 'CutoffHighPass_xFWHM: '+strcompress(-99,/rem)
		  endif
		  if (filterflag[i] eq '1') then begin
		    printf, lun, 'Spatial Filter: '+filtmethod[i]
		    printf, lun, 'HighPass [px]: '+sigfig(highpass,3)
		    printf, lun, 'LowPass [px]: '+sigfig(lowpass,3)
		  endif
		  if (method eq 'v' or method eq 'vc') then printf, lun, 'KLIP: variable'
		  if (method eq 'p' or method eq 'c') then printf, lun, 'KLIP: '+strcompress(klip[j],/rem)
		  if (method eq 'loci') then begin
		    printf, lun, 'Truncation: '+strcompress(klip[j],/rem)
		    printf, lun, 'Na: '+sigfig(na,3)
		    printf, lun, 'g: '+sigfig(g,3)
		  endif
		  if (method eq 'l') then begin
		    printf, lun, 'Rank: '+strcompress(klip[j],/rem)
		    printf, lun, 'Threshold: '+strcompress(uint(thresh),/rem)
		    printf, lun, 'Max. Iter: '+strcompress(uint(max_iter),/rem)
		  endif
		  if (method eq 'sa') then printf, lun, 'DistanceDelta [arcsec]: '+sigfig(distdelta,2)
		  if (method eq 'c') then printf, lun, 'ZonePCA/px: '+strcompress(zone[k],/rem)
		  if (method eq 'p' or method eq 'loci') then printf, lun, 'Step/FWHM: '+sigfig(step_init[k],2)
		  printf, lun, 'R_in/FWHM: '+sigfig(rin_init[l],2)
		  if (method eq 'p' or method eq 'loci') then printf, lun, 'Delta: '+sigfig(delta[m],2)
		  if (method eq 'v' or method eq 'vc') then printf, lun, 'Delta: variable'
		  ;printf, lun, ''
		  ;printf, lun, 'RightHanded: '+strmid(sigfig(right_handed,1),0,1)

		close, lun
		free_lun, lun

		;========================================================================

		if (fcpflag eq '1') then begin

		  print, ''
		  print, 'Fake Companion Injection'
		  ;print, ''

		  ;========================================================================

		  ; Fake companion throughput calibration parameters
		  ; ------------------------------------------------

		  sep_fcp = 4. ; radial angular separation between fake companions in the radial direction [in fwhm] 
		  n_rad = ROUND(dim/2/(sep_fcp*fwhm)) > 4  ;  number of fake companions along the radius (minimum 4)
		  n_br = 3 ; Number of branches of the pattern used to place our fake companions
		  rad_c_in = (rin_init[l]+0.5 > 2.0)*fwhm*plsc ; put the first FCP at 0.5 lambda/D from the PCA inner rim, but in any case beyond 2 lambda/D
		  rad_c_out_min = dim/2*plsc - 4.0*fwhm*plsc ; outer FCP at 4 lambda/D from the edge to prevent side effects
		  rad_c = rad_c_in + findgen(n_rad)*(rad_c_out_min-rad_c_in)/(n_rad-1)
		  pseudo_sigma_fcp = 30 ; SNR level at which the fake companions will be introduced

		  ;for LOCI another definition is needed because correction is not over full dim
		  if (method eq 'loci') then begin

		    sep_fcp = 4. 
		    n_br = 3
		    rad_c_in = (rin_init[l]+0.5 > 2.0)*fwhm*plsc
		    dr = floor (sqrt(!dpi*g*Na*fwhm^2/4)); radial "thickness" of the optimization region
		    rin = rin_init[l]*fwhm ; inner radius of the LOCI
		    rout = dim/2-1 - dr ; outer radius max (limited by box O size)
		    n_rad = ROUND(rout/(sep_fcp*fwhm)) > 4
		    step_law = 1.3
		    nstep = floor((rout-rin)^(1d/step_law)/(step_init[k]*fwhm))
		    rad_c_out_min = rin*plsc + floor((nstep*step_init[k]*fwhm)^step_law)*plsc- 4.0*fwhm*plsc	; outer FCP at 4 lambda/D from the edge to prevent side effects
		    rad_c = rad_c_in + findgen(n_rad)*(rad_c_out_min-rad_c_in)/(n_rad-1)
		    pseudo_sigma_fcp = 30 ; SNR level at which the fake companions will be introduced

		  endif
		  ;========================================================================

		  ;Fake companion injection at sigma * residual noise

		  if (recenter eq 'y') then tg_name_dyn = tg_name_ori+'_rsz_scaled_shift' else tg_name_dyn = tg_name_ori+'_rsz_scaled'

		  ;Injection of fake companions
		  ;PSF is currently not filtered! Disabled in fcp_am.pro
		  psffilt = 0
		  fcp_am, fwhm, transm_sci[xx], transm_fl[xx], ditfact, pseudo_sigma_fcp, rad_c, n_br, plsc, dim, camera_filter[xx], right_handed, filtmethod[i], PSFFILT=psffilt, /ADI, DISPLAY=display
    ; 	      if (filterflag[i] eq '1') then begin
    ; 
    ; 		psffilt = 1
    ; 		fcp_am, fwhm, transm_sci[xx], transm_fl[xx], ditfact, pseudo_sigma_fcp, rad_c, n_br, plsc, dim, camera_filter[xx], right_handed, filtmethod[i], PSFFILT=psffilt, /ADI, DISPLAY=display
    ; 
    ; 	      endif else begin
    ; 
    ; 		psffilt = 0
    ; 		fcp_am, fwhm, transm_sci[xx], transm_fl[xx], ditfact, pseudo_sigma_fcp, rad_c, n_br, plsc, dim, camera_filter[xx], right_handed, filtmethod[i], PSFFILT=psffilt, /ADI, DISPLAY=display
    ; 
    ; 	      endelse

		  print, ''

		  ;Spatial Filter
		  if (filterflag[i] eq '1') then MFILTER_V4b, filtmethod[i], HIGHPASS=highpass, LOWPASS=lowpass, DISPLAY=display;, /CA;/gauss;, /CA;, /MEDIAN_F

		  if (method eq 'c') then reduced_fcp_img_median = PCA_ADI_V39_am(cx, cy, dim, dim_init, fwhm, klip[j], zone[k], rin_init[l], delta[m], right_handed, 0, 0, NORMAL=normal, DISPLAY=display, /verbose)
		  if (method eq 'v') then reduced_fcp_img_median = SPCA_ADI_VAR(cx, cy, dim, dim_init, fwhm, plsc, klip[j], rin_init[l], step_init[k], delta[m], right_handed, 0, 0, pointdisk, 0, 0, /verbose, section=section)
		  if (method eq 'vc') then reduced_fcp_img_median = SPCA_ADI_VAR_CORR(cx, cy, dim, dim_init, fwhm, plsc, klip[j], rin_init[l], step_init[k], delta[m], right_handed, 0, 0, pointdisk, 0, 0, /verbose, /corr_file, section=section)
		  if (method eq 'p' and smethod eq '3') then reduced_fcp_img_median = SPCA_ADI_V13(cx, cy, dim, dim_init, fwhm, plsc, klip[j], rin_init[l], step_init[k], delta[m], right_handed, 0, 0, NORMAL=normal, /VERBOSE, DISPLAY=display)
		  if (method eq 'p' and smethod eq '5') then reduced_fcp_img_median = SPCA_ADI_V15(cx, cy, dim, dim_init, fwhm, plsc, klip[j], rin_init[l], step_init[k], delta[m], right_handed, 0, 0, NORMAL=normal, /VERBOSE, DISPLAY=display)
		  if (method eq 'l') then reduced_fcp_img_median = llsg_am(cx, cy, dim, dim_init, fwhm, plsc, right_handed, klip[j], thresh, max_iter, low_rank_mode, thresh_mode, rin_init[l], delta[m], 0, 0)
		  if (method eq 'loci') then reduced_fcp_img_median = loci_adi_v19(cx, cy, dim, dim_init, plsc, fwhm, rin_init[l], step_init[k], klip[j], Na, g, delta[m], right_handed, 0, 0, NORMAL=normal, /verbose, /svd);, DISPLAY=display
		  ;/svd, /nnls, /bvls, /damped)
		  if (method eq 'a') then reduced_fcp_img_median = adi_v01(cx, cy, dim, dim_init, plsc, fwhm, rin_init[l], right_handed, recenter, 0, 0, /verbose)
		  if (method eq 'sa') then reduced_fcp_img_median = sadi(cx, cy, dim, dim_init, plsc, fwhm, rin_init[l], delta[m], distdelta, right_handed, recenter, 0, 0, distdelta, /verbose)

		  ;Analyze fake companion throughput => alphar(r)
		  alpha_r_am, fwhm, reduced_img_median, n_br, rad_c, dim

		  ;Normalisation of the contrast curve : 5N(r)/alpha(r)
		  if keyword_set(sphere) then noise_curve_am, 0, dim/2, dim/2, fwhm, plsc, dim, sigma, rin_init[l],  LEVEL_FILE=procdir+'vec_'+tg_name_ori+'_level_r.fits', /ALPHA;, box=1;, /GAUSSFILT
		  if keyword_set(naco) then noise_curve_am, 0, dim/2, dim/2, fwhm, plsc, dim, sigma, rin_init[l],  LEVEL_FILE=procdir+'vec_'+tg_name_ori+'_level_r.fits', /ALPHA;, box=1;, /GAUSSFILT, agpm=agpm
		  if keyword_set(nici) then noise_curve_am, 0, dim/2, dim/2, fwhm, plsc, dim, sigma, rin_init[l],  LEVEL_FILE=procdir+'vec_'+tg_name_ori+'_level_r.fits', /ALPHA
		  if keyword_set(lmir) then noise_curve_am, 0, dim/2, dim/2, fwhm, plsc, dim, sigma, rin_init[l],  LEVEL_FILE=procdir+'vec_'+tg_name_ori+'_level_r.fits', /ALPHA
		  if keyword_set(nirc2) then noise_curve_am, 0, dim/2, dim/2, fwhm, plsc, dim, sigma, rin_init[l],  LEVEL_FILE=procdir+'vec_'+tg_name_ori+'_level_r.fits', /ALPHA

		  ;compute contrast
		  contrast_curve_am, fwhm, n_br, dim, rin_init[l], plsc, rad_c, sigma

    ; 	      if keyword_set(pmodel) then begin; (n_elements(stellmag) gt 0 and n_elements(age) gt 0) then begin
    ; 
    ; 		if (strmatch(camera_filter[xx], '*J*') eq 1) then mag = (stellmag[0])[0]
    ; 		if (strmatch(camera_filter[xx], '*H*') eq 1) then mag = (stellmag[1])[0]
    ; 		if (strmatch(camera_filter[xx], '*K*') eq 1) then mag = (stellmag[2])[0]
    ; 		if (strmatch(camera_filter[xx], 'Lp') eq 1) then mag = (stellmag[0])[0]
    ; 
    ; 		;compute mass detection limit based on BTsettl models
    ; 		print, ''
    ; 		print, 'Estimating Mass Detection Limits...'
    ; 
    ; 		if (finite(plx) ne 1) then plx = 1./100.
    ; 
    ; 		planet_lower_mass_limit, dim, plsc, sigma, mag, plx, age, camera_filter[xx], sphere=sphere, naco=naco
    ; 
    ; 	      endif

		  ;========================================================================

; 		  ;this is a check if we can recover e.g. 5sigma planets
; 		  ;at this point we have noise level, flux loss
; 
; 		  if (recenter eq 'y') then tg_name_dyn = tg_name_ori+'_rsz_scaled_shift' else tg_name_dyn = tg_name_ori+'_rsz_scaled'
; 
; 		  ;PSF is currently not filtered! Disabled in fcp_am.pro
; 		  fcp_am, fwhm, transm_sci[xx], transm_fl[xx], ditfact, sigma, rad_c, n_br, plsc, dim, camera_filter[xx], right_handed, filtmethod[i], PSFFILT=psffilt, /ADI, /CHECK;, DISPLAY=display;
; 
; 		  print, ''
; 
; 		  ;Spatial Filter
; 		  if (filterflag[i] eq '1') then MFILTER_V4b, filtmethod[i], HIGHPASS=highpass, LOWPASS=lowpass, DISPLAY=display;, /CA;/gauss;, /CA;, /MEDIAN_F
; 
; 		  ;if (method eq 'c') then reduced_fcp_img_median = PCA_ADI_V39(cx, cy, dim, dim_init, fwhm, klip[j], zone[k], rin_init[l], right_handed, NORMAL=normal, DISPLAY=display, /verbose)
; 		  if (method eq 'v') then reduced_fcp_img_median = SPCA_ADI_VAR(cx, cy, dim, dim_init, fwhm, plsc, klip[j], rin_init[l], step_init[k], delta[m], right_handed, 0, 0, pointdisk, 0, 0, /verbose, section=section)
; 		  if (method eq 'p' and smethod eq '3') then reduced_fcp_img_median = SPCA_ADI_V13(cx, cy, dim, dim_init, fwhm, plsc, klip[j], rin_init[l], step_init[k], delta[m], right_handed, 0, 0, NORMAL=normal, /VERBOSE, DISPLAY=display)
; 		  if (method eq 'p' and smethod eq '5') then reduced_fcp_img_median = SPCA_ADI_V15(cx, cy, dim, dim_init, fwhm, plsc, klip[j], rin_init[l], step_init[k], delta[m], right_handed, 0, 0, NORMAL=normal, /VERBOSE, DISPLAY=display)
; 		  if (method eq 'l') then reduced_fcp_img_median = llsg_am(cx, cy, dim, dim_init, fwhm, plsc, right_handed, klip[j], thresh, max_iter, low_rank_mode, thresh_mode, rin_init[l], delta[m], 0, 0)
; 		  if (method eq 'loci') then reduced_fcp_img_median = loci_adi_v19(cx, cy, dim, dim_init, plsc, fwhm, rin_init[l], step_init[k], klip[j], Na, g, delta[m], right_handed, 0, 0, NORMAL=normal, /verbose, /svd);, DISPLAY=display
; 		  ;/svd, /nnls, /bvls, /damped)
; 		  if (method eq 'a') then reduced_fcp_img_median = adi_v01(cx, cy, dim, dim_init, plsc, fwhm, rin_init[l], right_handed, recenter, 0, 0)

		  ;========================================================================

		  ; Plot the 2D contrast map, based on a 3-pixel box
		  alpha_2d = READFITS(procdir+'img_'+tg_name_bas+'_alpha_r_2D.fits',header, /silent)
		  alpha_2d[where(alpha_2d eq 0)] = 1D0
		  noise_map = image_stddev(filter_image(reduced_img_median, fwhm_gaussian=0.5*fwhm), 1.5*fwhm)/(alpha_2d>0.01)  ; use half the box size

		  idx = where(reduced_img_median eq 0.d0)
		  if (idx[0] ne -1) then begin

		    idx0 = array_indices(reduced_img_median, idx)
		    if (idx0[0] ne -1) then begin

		      for yy=0L,n_elements(idx0[0,*])-1 do noise_map[idx0[0,yy],idx0[1,yy]] = !values.d_nan

		    endif

		  endif

		  ;--------------------------------------------------------------------------

		  mask_t = shift(dist(dim), dim/2., dim/2.)
		  mask = mask_t ge fwhm*rin_init[l]

		  ;--------------------------------------------------------------------------
		  ;Contrast Map

		  psf = readfits(procdir+'img_'+tg_name_bas+'_scaledPSF.fits', /silent)
		  maxpsf = (max(psf))[0]

		  cmap = alog10(sigma*noise_map/maxpsf)*mask

		  set_plot, 'ps'
		  fn = datadir+'plot_'+tg_name_bas+'_contrast_map.ps'
		  device, filename=fn, /color, XSIZE=35, YSIZE=30
		  !p.font=0
		  !p.thick=3
		  !x.thick=3
		  !y.thick=3
		  loadct, 70, /silent

		  thisLetter = "104B
		  deltap = '!9' + String(thisLetter) + '!X'
		  thisLetter = "163B
		  sigmal = '!9' + String(thisLetter) + '!X'
		  cgimage, cmap, /scale, xr=[dim/2,-dim/2]*plsc, yr=[-dim/2,dim/2]*plsc, $
			xtitle=deltap+'RA [arcsec]', ytitle=deltap+'DEC [arcsec]', /axis, charsize=2, position=[0.09,0.08,0.85,0.98], minvalue=min(cmap,/nan), maxvalue=max(cmap,/nan) ;, background=cgcolor('white')
		  oplot, [-dim/2,dim/2]*plsc, [0,0]*plsc, linestyle=1 & oplot, [0,0]*plsc, [-dim/2,dim/2]*plsc, linestyle=1
		  cgcolorbar, division=5, range=[min(cmap,/nan),max(cmap,/nan)], /vertical, /right, position=[0.86, 0.08, 0.9, 0.98], title=STRING(sigma, format='(I0)')+sigmal+' Contrast (log)', charsize=2, ticklen=1;, tickinterval=0.2

		  loadct, 0, /silent
		  !p.thick=1
		  !x.thick=1
		  !y.thick=1

		  device,/close
		  set_plot,'x'
		  spawn, '/home/amueller/work/IDLlibs/epstopdf.pl '+fn
		  spawn, 'rm '+fn

		  writefits, datadir+tg_name_bas+'_contrast_map.fits', cmap

		  ;--------------------------------------------------------------------------
		  ;SNR Map

		  if (filterflag[i] eq '0') then tmp = filter_image(reduced_img_median, fwhm_gaussian=0.5*fwhm) else tmp = reduced_img_median
		  smap = tmp/rmsmap_am(tmp, fwhm)*mask
		  ;smap = filter_image(reduced_img_median, fwhm_gaussian=fwhm)/noise_map
		  set_plot, 'ps'
		  fn = datadir+'plot_'+tg_name_bas+'_SNR_map.ps'
		  device, filename=fn, /color, XSIZE=35, YSIZE=30
		  !p.font=0
		  !p.thick=3
		  !x.thick=3
		  !y.thick=3
		  loadct, 70, /silent

		  thisLetter = "104B
		  deltap = '!9' + String(thisLetter) + '!X'
		  cgimage, smap, /scale, xr=[dim/2,-dim/2]*plsc, yr=[-dim/2,dim/2]*plsc, $
			xtitle=deltap+'RA [arcsec]', ytitle=deltap+'DEC [arcsec]', /axis, charsize=2, position=[0.09,0.08,0.85,0.98], minvalue=-5, maxvalue=5
		  oplot, [-dim/2,dim/2]*plsc, [0,0]*plsc, linestyle=1 & oplot, [0,0]*plsc, [-dim/2,dim/2]*plsc, linestyle=1
		  cgcolorbar, division=3, range=[-5,5], /vertical, /right, position=[0.86, 0.08, 0.9, 0.98], title='SNR', charsize=2, tickinterval=1, ticklen=1

		  loadct, 0, /silent
		  !p.thick=1
		  !x.thick=1
		  !y.thick=1

		  device,/close
		  set_plot,'x'
		  spawn, '/home/amueller/work/IDLlibs/epstopdf.pl '+fn
		  spawn, 'rm '+fn

		  writefits, datadir+tg_name_bas+'_SNR_map.fits', smap

		  ;--------------------------------------------------------------------------
		  ;SNR FCP Map

		  if (filterflag[i] eq '0') then tmp = filter_image(reduced_fcp_img_median, fwhm_gaussian=0.5*fwhm) else tmp = reduced_fcp_img_median
		  sfcpmap = tmp/rmsmap_am(tmp, fwhm)*mask
    ; 	      sfcpmap = filter_image(reduced_fcp_img_median, fwhm_gaussian=fwhm)/noise_map
		  set_plot, 'ps'
		  fn = datadir+'plot_'+tg_name_bas+'_SNR_fcp_map.ps'
		  device, filename=fn, /color, XSIZE=35, YSIZE=30
		  !p.font=0
		  !p.thick=3
		  !x.thick=3
		  !y.thick=3
		  loadct, 70, /silent

		  thisLetter = "104B
		  deltap = '!9' + String(thisLetter) + '!X'
		  cgimage, sfcpmap, /scale, xr=[dim/2,-dim/2]*plsc, yr=[-dim/2,dim/2]*plsc, $
			xtitle=deltap+'RA [arcsec]', ytitle=deltap+'DEC [arcsec]', /axis, charsize=2, position=[0.09,0.08,0.85,0.98], minvalue=-5, maxvalue=5
		  oplot, [-dim/2,dim/2]*plsc, [0,0]*plsc, linestyle=1 & oplot, [0,0]*plsc, [-dim/2,dim/2]*plsc, linestyle=1
		  cgcolorbar, division=3, range=[-5,5], /vertical, /right, position=[0.86, 0.08, 0.9, 0.98], title='SNR', charsize=2, tickinterval=1, ticklen=1

		  loadct, 0, /silent
		  !p.thick=1
		  !x.thick=1
		  !y.thick=1

		  device,/close
		  set_plot,'x'
		  spawn, '/home/amueller/work/IDLlibs/epstopdf.pl '+fn
  		  spawn, 'rm '+fn

		  writefits, datadir+tg_name_bas+'_SNR_fcp_map.fits', sfcpmap

		  ;========================================================================

		endif

		;========================================================================

		;move results to resultsdir

		spawn, 'cp '+datadir+filepsf[xx]+' '+resdir
		spawn, 'rm '+datadir+'*_rsz_dc.fits'	;remove the resized image to save some space
		spawn, 'mv '+datadir+'*_rsz_* '+resdir
		if (fcpflag eq '1') then spawn, 'mv '+datadir+'*_level_*.fits '+resdir
		spawn, 'mv '+datadir+'setup*.txt '+resdir
		if (canddet eq '1') then spawn, 'mv '+datadir+star+'*.sav '+resdir
                spawn, 'mv '+datadir+'*_scaledPSF*.fits '+resdir
		if (method eq 'v' or method eq 'vc') then spawn, 'mv '+datadir+'Annulus_*.txt '+resdirkd
		if (method eq 'vc') then spawn, 'mv '+datadir+'*.sav '+resdir

; 		if (method eq 'p') then begin
; 
; 		  fev = file_search(datadir+'PCA_Eigenvalues*.txt', count=nev)
; 		  for xxx=0,nev-1 do spawn, 'mv '+fev[xxx]+' '+resdirev
; 
; 		endif

		if (fcpflag eq '1') then begin

		  spawn, 'mv '+datadir+'*_map.* '+resdir
		  ;spawn, 'mv '+datadir+'*.ps '+resdir
		  ;spawn, 'mv '+datadir+'*.pdf '+resdir
		  spawn, 'mv '+datadir+'*alpha*.fits '+resdir
		  spawn, 'mv '+datadir+'*ratio*.fits '+resdir
		  ;spawn, 'mv '+datadir+'*aperture*.fits '+resdir
		  spawn, 'mv '+datadir+'*_coords.fits '+resdir
		  spawn, 'mv '+datadir+'*_fcparray.fits '+resdir
		  spawn, 'mv '+datadir+'*_contrast* '+resdir
		  ;spawn, 'mv '+datadir+'*_phot* '+resdir

		endif

	      endif

	      ;========================================================================



	      ;create summary

	      if (fcpflag eq '1' and selins ne 'l') then begin

		;Candidates detected
		if (canddet eq '1') then restore, resdir+star+'_'+camera_filter[xx]+'_Candidates.sav'


		if (filterflag eq '1') then imfile = file_search(resdir+'*'+camera_filter[xx]+'*scaled*median.fits') else imfile = file_search(resdir+'*'+camera_filter[xx]+'*scaled*median_convolved.fits')
		idx1 = strmatch(imfile, '*_fcp_*')
		idx2 = where(idx1 eq 0)
		if (idx2[0] eq -1) then begin

		  print, 'No final image found. Stop.'
		  stop

		endif

		if (n_elements(idx2) gt 1) then begin

		  print, 'No unique file found. Stop.'
		  stop

		endif

		imfile = (imfile[idx2])[0]
		im = mrdfits(imfile, 0, /silent)

		if (strmatch(camera_filter[xx], '*J*') eq 1) then mag = (stellmag[0])[0]
		if (strmatch(camera_filter[xx], '*PaB*') eq 1) then mag = (stellmag[0])[0]  ;use J mag for the moment
		if (strmatch(camera_filter[xx], '*H*') eq 1) then mag = (stellmag[1])[0]
		if (strmatch(camera_filter[xx], '*K*') eq 1) then mag = (stellmag[2])[0]
		if (strmatch(camera_filter[xx], 'Lp') eq 1) then mag = (stellmag[0])[0]
		if (strmatch(camera_filter[xx], 'Ks') eq 1) then mag = (stellmag[1])[0]
		if (strmatch(camera_filter[xx], 'dx') eq 1) then mag = (stellmag[0])[0]
		if (strmatch(camera_filter[xx], 'sx') eq 1) then mag = (stellmag[0])[0]
		
		contrast = mrdfits(resdir+'vec_'+camera_filter[xx]+'_contrast_r.fits', 0, /silent)

		detection_mag = logm(1./contrast)+mag
		detection = 1./exp(detection_mag*alog(10.d0^0.4d0))
		dist = dindgen(n_elements(contrast))*plsc

		;find RAW data to get header info
		if (selins eq 'n') then begin

		  tmpfile = file_search(path+'../Reduced/'+'cube_*reduced.fits', count=nall)
		  hdrtmp = headfits(tmpfile[0], exten=0, /silent)
		  mjdtmp = double(get_eso_keyword(hdrtmp, 'MJD-OBS'))

		  if (mjdtmp eq 0.d0) then scifile = file_search(path+'NACO*fits', count=nall) else scifile = tmpfile

		endif
		;if (selins eq 'n') then scifile = file_search(path+'NACO*fits', count=nall)

		if (selins eq 'ni') then begin

		  tmpfile = file_search(path+'../Reduced/'+camera_filter[xx]+'_cube_*reduced.fits', count=nall)
		  hdrtmp = headfits(tmpfile[0], exten=0, /silent)
		  mjdtmp = double(get_eso_keyword(hdrtmp, 'MJD_OBS'))
		  scifile = tmpfile

		endif

		if (selins eq 's') then begin

		  scifile = file_search(path+'SPH*fits', count=nall)

		endif

		catgall = strarr(nall)
		typeall = strarr(nall)
		for sc=0,nall-1 do begin

		  hdr = headfits(scifile[sc], exten=0, /silent)
		  catgall[sc] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DPR CATG'),/rem)
		  typeall[sc] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DPR TYPE'),/rem)

		endfor

		if (selins eq 'n') then idx = where(catgall eq 'SCIENCE')
  ; 	      if (selins eq 's') then idx = where(catgall eq 'SCIENCE' and typeall eq 'OBJECT,CENTER')
		if (selins eq 's') then idx = where(catgall eq 'SCIENCE' and typeall eq 'OBJECT')

		scifile = scifile[idx]
		nsci = n_elements(scifile)

		seeing = dblarr(nsci)
		alt = dblarr(nsci)
		rh = dblarr(nsci)
		tau = dblarr(nsci)
		wind = dblarr(nsci)
		lst  = dblarr(nsci)
		ha = dblarr(nsci)
		mjd = dblarr(nsci)
		pupilpos = dblarr(nsci)
		;dit = dblarr(nsci)

		for sc=0,nsci-1 do begin

		  hdr = headfits(scifile[sc], exten=0, /silent)

		  if (selins eq 'n' or selins eq 's') then begin

		    if (sc eq 0) then progid = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO OBS PROG ID'),/rem)
		    if (sc eq 0) then ra = double(get_eso_keyword(hdr, 'RA'))/15.
		    lst[sc] = double(get_eso_keyword(hdr, 'LST'))/3600.
		    ha[sc] = lst[sc] - ra
		    mjd[sc] = double(get_eso_keyword(hdr, 'MJD-OBS'))+2400000.d0	;without 0.5 on purpose to get the night right
		    alt[sc] = double(get_eso_keyword(hdr, 'HIERARCH ESO TEL ALT'))
		    rh[sc] = double(get_eso_keyword(hdr, 'HIERARCH ESO TEL AMBI RHUM'))
		    if (selins eq 'n') then tau[sc] = double(get_eso_keyword(hdr, 'HIERARCH ESO TEL AMBI TAU0'))*1.d3
		    if (selins eq 's') then tau[sc] = double(get_eso_keyword(hdr, 'HIERARCH ESO TEL AMBI TAU0'))*1.d3
		    wind[sc] = double(get_eso_keyword(hdr, 'HIERARCH ESO TEL AMBI WINDSP'))
		    seeing[sc] = double(get_eso_keyword(hdr, 'HIERARCH ESO TEL IA FWHM'))
		    if (sc eq 0) then lat = double(get_eso_keyword(hdr, 'HIERARCH ESO TEL GEOLAT'))
		    ;HIERARCH ESO TEL IA FWHMLIN
		    ;HIERARCH ESO TEL IA FWHMLINOBS
		    ;dit[sc] = double(get_eso_keyword(hdr, 'HIERARCH ESO DET DIT'))

		    if (selins eq 'n') then pupilpos[sc] = double(get_eso_keyword(hdr, 'HIERARCH ESO ADA PUPILPOS'))

		  endif

		  if (selins eq 'ni') then begin

		    if (sc eq 0) then progid = strcompress(get_eso_keyword(hdr, 'GEMPRGID'),/rem)
		    if (sc eq 0) then ra = double(get_eso_keyword(hdr, 'RA'))/15.
		    tmp = get_eso_keyword(hdr, 'ST')
		    tmp1 = double(strmid(tmp,0,2))
		    tmp2 = double(strmid(tmp,3,2))
		    tmp3 = double(strmid(tmp,6,4))
		    lst[sc] = tmp1+tmp2/60.+tmp3/3600.

		    ha[sc] = lst[sc] - ra
		    mjd[sc] = double(get_eso_keyword(hdr, 'MJD_OBS'))+2400000.d0	;without 0.5 on purpose to get the night right
		    alt[sc] = double(get_eso_keyword(hdr, 'ELEVATIO'))
		    rh[sc] = double(get_eso_keyword(hdr, 'HUMIDITY'))
		    tau[sc] = -99.
		    wind[sc] = double(get_eso_keyword(hdr, 'WINDSPEE'))
		    seeing[sc] = -99.
		    lat = ten(-30.d0, 14.d0, 26.70d0)
		    iaa = double(get_eso_keyword(hdr, 'IAA'))
		    crpa = double(get_eso_keyword(hdr, 'CRPA'))

		  endif

		  if (selins eq 'n2') then begin

		    if (sc eq 0) then progid = strcompress(get_eso_keyword(hdr, 'PROGID'),/rem)
		    if (sc eq 0) then ra = double(get_eso_keyword(hdr, 'RA'))/15.
		    tmp = get_eso_keyword(hdr, 'ST')
		    tmp1 = double(strmid(tmp,0,2))
		    tmp2 = double(strmid(tmp,3,2))
		    tmp3 = double(strmid(tmp,6,5))
		    lst[sc] = tmp1+tmp2/60.+tmp3/3600.

		    ha[sc] = lst[sc] - ra
		    mjd[sc] = double(get_eso_keyword(hdr, 'MJD-OBS'))+2400000.d0	;without 0.5 on purpose to get the night right
		    alt[sc] = double(get_eso_keyword(hdr, 'EL'))
		    rh[sc] = double(get_eso_keyword(hdr, 'WXOUTHUM'))
		    tau[sc] = -99.
		    wind[sc] = double(get_eso_keyword(hdr, 'WXWNDSP'))
		    seeing[sc] = -99.
		    lat = ten(-30.d0, 14.d0, 26.70d0)
; 		    iaa = double(get_eso_keyword(hdr, 'IAA'))
; 		    crpa = double(get_eso_keyword(hdr, 'CRPA'))

		  endif
		  
		endfor

		tobs = (max(mjd)-min(mjd))*24.	;hours
		;dit_sci
		;dit_fl
		;mag

		paral = mrdfits(datadir+'vec_'+tg_name_bas+'_paral.fits', 0, /silent)
		haused = mrdfits(datadir+'vec_'+tg_name_bas+'_paral.fits', 1, /silent)

		if (n_elements(haused) gt 1) then begin

		  ;ha = haused
		  flagnewred2 = 1

		endif else begin

		  flagnewred2 = 0

		endelse

		;fieldrot = abs(min(paral)-max(paral))

		caldat, min(mjd), mm, dd, yy, t1, t2, t3
		if (dd lt 10.) then dd = '0'+strcompress(dd,/rem) else dd = strcompress(dd,/rem)
		if (mm lt 10.) then mm = '0'+strcompress(mm,/rem) else mm = strcompress(mm,/rem)
		yy = strcompress(yy,/rem)
		night = dd+'.'+mm+'.'+yy

		if keyword_set(naco) then begin

		  if (agpm eq 1) then begin

		    if (strmatch(resdir, '*Archive*') eq 1) then fn = resdir+'Archive_'+star+'_'+camera_filter[xx]+'_'+yy+mm+dd+'_overview'+dirflag $
		      else fn = resdir+star+'_'+camera_filter[xx]+'_'+yy+mm+dd+'_overview'+dirflag

		  endif else begin

		    if (strmatch(resdir, '*Archive*') eq 1) then fn = resdir+'Archive_'+star+'_noAGPM_'+camera_filter[xx]+'_'+yy+mm+dd+'_overview'+dirflag $
		      else fn = resdir+star+'_noAGPM_'+camera_filter[xx]+'_'+yy+mm+dd+'_overview'+dirflag

		  endelse

		endif

		if keyword_set(nici) then fn = resdir+star+'_'+camera_filter[xx]+'_'+yy+mm+dd+'_overview'+dirflag

		if keyword_set(sphere) then fn = resdir+star+'_'+camera_filter[xx]+'_'+yy+mm+dd+'_overview'+dirflag

		;---------------------------------------------------------------------------------

		loadct, 1, /silent

		maxexp = 10.^(double(ceil(abs(alog10(max(detection, /nan)))))-1.)
		minexp = 10.^double(ceil(abs(alog10(min(detection, /nan)))))

		!p.font=0
		!p.thick=3
		!x.thick=3
		!y.thick=3
		;!p.multi=[0,1,0]

		thisLetter = "104B
		deltap = '!9' + String(thisLetter) + '!X'
		thisLetter = "163B
		sigmal = '!9' + String(thisLetter) + '!X'
		thisLetter = "164B
		taul = '!9' + String(thisLetter) + '!X'
		thisLetter = "142B
		betal = '!9' + String(thisLetter) + '!X'
		quote = "'"

		set_plot, 'ps'
		device, filename=fn+'.ps', /color, XSIZE=21, YSIZE=29.7

; 		  cgimage, im, /axis, xr=[dim/2,-dim/2]*plsc, yr=[-dim/2,dim/2]*plsc, xtitle=deltap+'RA [arcsec]', ytitle=deltap+'DEC [arcsec]', position=[0.11,0.48,0.8,0.96], charsize=1., stretch='clip', clip=0.1;, minvalue=-5.d-4, maxvalue=5.d-4;;, minvalue=-10, maxvalue=10;, minvalue=min(im), maxvalue=max(im), title=star
		  mask_t = shift(dist(dim), dim/2., dim/2.)
		  mask = mask_t ge rin_init[l]*fwhm
		  im = im*mask
		  cgimage, scale_image_am(im), /axis, xr=[dim/2,-dim/2]*plsc, yr=[-dim/2,dim/2]*plsc, xtitle=deltap+'RA [arcsec]', ytitle=deltap+'DEC [arcsec]', position=[0.11,0.48,0.8,0.96], charsize=1.
		  oplot, [-dim/2,dim/2]*plsc, [0,0]*plsc, linestyle=1, color=cgcolor('gray')
		  oplot, [0,0]*plsc, [-dim/2,dim/2]*plsc, linestyle=1, color=cgcolor('gray')

		  if (canddet eq '1') then begin

		    if (ncand ne 0) then begin


		      for cc=0,ncand-1 do begin

			x0 = ((dim/2.-xf[cc])-2.*fwhm)*plsc
			x1 = ((dim/2.-xf[cc])+2.*fwhm)*plsc
			y0 = ((yf[cc]-dim/2.)-2.*fwhm)*plsc
			y1 = ((yf[cc]-dim/2.)+2.*fwhm)*plsc

			plots, [x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], /data, color=cgcolor('green'), thick=4

			xyouts, x1+(0.1*abs(x0-x1)), y1, strcompress(cc+1), /data, alignment=1, color=cgcolor('green'), charsize=0.7

		      endfor

		      ;for cc=0,ncand-1 do tvcircle, /data, 2*fwhm, (dim/2.-xf[cc])*plsc, (yf[cc]-dim/2.)*plsc, color=cgcolor('green'), thick=2

		    endif

		  endif


		  cgplot, /noerase, dist, logm(1./contrast), xtitle='Angular separation [arcsec]', ytitle=strcompress(uint(sigma),/rem)+sigmal+' Contrast ['+deltap+'mag]', xr=[0., ceil(max(dist)*10.)/10.], xst=1, yst=1, XTicklen=1.0, YTicklen=1.0, XGridStyle=1, YGridStyle=1, yr=[ceil(max(logm(1./contrast), /nan)),floor(min(logm(1./contrast), /nan))-0.5], position=[0.11,0.28,0.8,0.43], charsize=1., color=cgcolor('blue');, /nodata
    ; 	      oplot, dist, detection_mag, thick=5, color=cgcolor('blue')

		  ;idx = where(ha ne 0.)
		  ;ha = ha[idx]
		  ;alt = alt[idx]
		  idxalt = where(alt ne 0.)
		  ha = ha[idxalt]
		  alt = alt[idxalt]
		  

		  idxha = where(ha lt -12.)
		  if (idxha[0] ne -1) then ha[idxha] = ha[idxha]+24.
		  idxha = where(ha gt 12.)
		  if (idxha[0] ne -1) then ha[idxha] = ha[idxha]-24.


  ; 		plot, ha, alt, position=[0.35,0.05,0.8,0.22], /normal, /noerase, charsize=0.7, xtitle='Hour Angle [hr]', ytitle='Altitude [deg]', /yn, psym=sym(1), symsize=0.5
		  cgplot, ha, alt, position=[0.35,0.05,0.8,0.22], /normal, /noerase, charsize=0.7, xtitle='Hour Angle [hr]', ytitle='Altitude [deg]', /yn, psym=sym(1), symsize=0.5, yst=8
		  cgAxis, YAxis=1, /Save, charsize=0.7, ytitle='Parallactic Angle [deg]', color=cgcolor('blue'), yst=1, yr=[-180,180];yr=[floor(min(parang)), ceil(max(parang))]

		  ;plotting the parallactic angle but only for frames used! Therefore, the grapg gives a good overview what we have for data in total (black points for altitude).
		  decdec = ten(double(strmid(dec,0,3)), double(strmid(dec,3,2)), double(strmid(dec,5,6)))

		  if (flagnewred2 eq 0) then begin

		    if (selins eq 'n') then begin

		      haparang = linspace(!x.crange[0], !x.crange[1], 200)
		      parang = parangle(haparang, decdec, lat)

		      if (decdec gt lat) then begin

			for cc=0,n_elements(haparang)-1 do begin

			  if (parang[cc] lt 0.) then parang[cc] = parang[cc]+360.d0

			endfor

		      endif

		    endif else begin

		      parang = parangle(ha, decdec, lat)

		      if (decdec gt lat) then begin

			for cc=0,n_elements(haparang)-1 do begin

			  if (parang[cc] lt 0.) then parang[cc] = parang[cc]+360.d0

			endfor

		      endif

		    endelse


		    if (selins eq 'n') then begin

		      ;print, 'Numbner of Pupilpos and parang are not the same! Stop.'
		      if (abs(min(pupilpos)-max(pupilpos)) gt 2.) then begin

			print, 'change of Pupil position! Stop.'
			stop

		      endif

		      parang = parang-(90.d0+(89.44-pupilpos[0]))
		      ;I use only one pupil position because it is read from the reduced frames. But paral can vary becasue of frame selection. If pupilpos is more or less constant its fine. In the case of P97/201606/HD163296 pupilpos changed significantly because of a wrong setup

		    endif

		    if (selins eq 'n2') then parang = parang
		    if (selins eq 'ni') then parang = parang+iaa+crpa

		    if (selins eq 's') then parang = parang+1.75+135.99

  ; 		  cgoplot, haparang, parang, color=cgcolor('green');, psym=sym(1), symsize=0.5

		    if (selins eq 'n') then begin

		      plothaparal = dblarr(n_elements(paral))
		      for cc=0,n_elements(paral)-1 do begin

			idx = closest(parang, paral[cc])
			plothaparal[cc] = haparang[idx]

		      endfor

		    endif else begin

		      plothaparal = ha
		      paral = parang

		    endelse

		    fieldrot1 = abs(min(paral)-max(paral))

		    idx = where(paral gt 180.)
		    if (idx[0] ne -1) then paral[idx] = paral[idx]-360.
		    idx = where(paral lt -180.)
		    if (idx[0] ne -1) then paral[idx] = paral[idx]+360.

		    fieldrot2 = abs(min(paral)-max(paral))	;i have to do this shit because of the wrapping shit
		    if (fieldrot1 le fieldrot2) then fieldrot = fieldrot1 else fieldrot = fieldrot2

		  endif else begin

		    fieldrot1 = abs(min(paral)-max(paral))
		    plothaparal = haused
		    idx = where(paral gt 180.)
		    if (idx[0] ne -1) then paral[idx] = paral[idx]-360.
		    idx = where(paral lt -180.)
		    if (idx[0] ne -1) then paral[idx] = paral[idx]+360.

		    fieldrot2 = abs(min(paral)-max(paral))
		    if (fieldrot1 le fieldrot2) then fieldrot = fieldrot1 else fieldrot = fieldrot2

		  endelse

		  cgoplot, plothaparal, paral, color=cgcolor('blue'), psym=sym(1), symsize=0.5


; 		  cgoplot, haparang, paral, color=cgcolor('blue'), psym=sym(1), symsize=0.5

		  ;plot, ha, seeing, position=[0.3,0.05,0.8,0.2], /normal, /noerase, charsize=0.7, xtitle='Hour Angle [hr]', ytitle='Seeing [asec]'

		  if (canddet eq '1') then begin

		    xyouts, 0.81, 0.95, '#   Dist.  PA    SNR   Qual.', alignment=0., charsize=0.7, color=cgcolor('black'), /normal
		    if (ncand gt 0 and ncand lt 50) then begin

		      for cc=0,ncand-1 do begin

			outcc = strcompress(uint(cc+1),/rem)+'  '+sigfig(rhosf[cc],4)+'  '+sigfig(pasf[cc],3)+'  '+sigfig(snsf[cc],3)+'  '+sigfig(fqsf[cc],2)
			xyouts, 0.81, 0.95-(cc+1)*0.01, outcc, alignment=0., charsize=0.7, color=cgcolor('black'), /normal

		      endfor

		    endif

		    if (ncand ge 50) then begin

		      for cc=0,49 do begin

			outcc = strcompress(uint(cc+1),/rem)+'  '+sigfig(rhosf[cc],4)+'  '+sigfig(pasf[cc],3)+'  '+sigfig(snsf[cc],3)+'  '+sigfig(fqsf[cc],2)
			xyouts, 0.81, 0.95-(cc+1)*0.01, outcc, alignment=0., charsize=0.7, color=cgcolor('black'), /normal

		      endfor

			xyouts, 0.81, 0.95-(cc+1)*0.01, '...', alignment=0., charsize=0.7, color=cgcolor('black'), /normal
			cc = cc+1
			xyouts, 0.81, 0.95-(cc+1)*0.01, '...', alignment=0., charsize=0.7, color=cgcolor('black'), /normal
			cc = cc+1
			xyouts, 0.81, 0.95-(cc+1)*0.01, '...', alignment=0., charsize=0.7, color=cgcolor('black'), /normal

		    endif

		  endif

		  if (star eq 'betaPic') then dumstar = betal+'Pic' else dumstar = star
		  xyouts, 0.5, 0.97, dumstar, color=cgcolor('black'), charsize=1.5, alignment=0.5, /normal

		  if (method eq 'c') then dumalgo = 'PCA Box'
		  if (method eq 'v') then dumalgo = 'PCA Var'
		  if (method eq 'vc') then dumalgo = 'PCA Var Corr'
		  if (method eq 'p') then dumalgo = 'PCA'
		  if (method eq 'l') then dumalgo = 'LLSG'
		  if (method eq 'a') then dumalgo = 'ADI'
		  if (method eq 'sa') then dumalgo = 'sADI'

		  camfiltdum = camera_filter[xx]
		  if (camfiltdum eq 'Lp') then camfiltdum = 'L'+quote
		  if (camfiltdum eq 'J2') then camfiltdum = 'J'
		  if (camfiltdum eq 'J3') then camfiltdum = 'J'
		  if (camfiltdum eq 'H2') then camfiltdum = 'H'
		  if (camfiltdum eq 'H3') then camfiltdum = 'H'
		  if (camfiltdum eq 'K1') then camfiltdum = 'K'
		  if (camfiltdum eq 'K2') then camfiltdum = 'K'
		  if (camfiltdum eq 'Ms') then camfiltdum = 'Ms'

		  if (selins eq 'n') then begin

		    if (agpm eq 1) then dumagpm = 'AGPM' else dumagpm = 'sat. PSF'
		    legend, [ $
		    string('Night', format='(a-26)')+string(': ',format='(a-5)')+string(night,format='(a-12)'), $
		    string('Prog. ID', format='(a-24)')+string(': ',format='(a-5)')+string(progid,format='(a-13)'), $
		    string('Obs. mode', format='(a-21)')+string(': ',format='(a-5)')+string(dumagpm,format='(a-12)'), $
		    string('Algorithm', format='(a-24)')+string(': ',format='(a-5)')+string(dumalgo,format='(a-12)'), $
		    string('Field rotation [deg]', format='(a-20)')+string(': ',format='(a-5)')+string(sigfig(fieldrot,4),format='(a-12)'), $
		    string('Obs. time [hr]', format='(a-22)')+string(': ',format='(a-5)')+string(sigfig(tobs,3),format='(a-12)'), $
		    string('DIT science [ms]', format='(a-19)')+string(': ',format='(a-5)')+string(sigfig(dit_sci*1.d3,4),format='(a-12)'), $
		    string('DIT flux [ms]', format='(a-23)')+string(': ',format='(a-5)')+string(sigfig(dit_fl*1.d3,4),format='(a-12)'), $
		    string('Seeing [asec]', format='(a-21)')+string(': ',format='(a-5)')+string(sigfig(mean(seeing),2),format='(a-12)'), $
		    string(taul+'!D0!N [ms]', format='(a-34)')+string(': ',format='(a-5)')+string(sigfig(mean(tau),3),format='(a-12)'), $
		    string('Wind speed [m/s]', format='(a-18)')+string(': ',format='(a-5)')+string(sigfig(mean(wind),3),format='(a-12)'), $
		    string('RH [%]', format='(a-24)')+string(': ',format='(a-5)')+string(sigfig(mean(rh),3),format='(a-12)'), $
		    string(camfiltdum+' [mag]', format='(a-25)')+string(': ',format='(a-5)')+string(sigfig(mag,4),format='(a-12)')], $
		    box=0, margin=0, /bottom, /left, position=[0.04, 0.09], /normal, charsize=0.7

		  endif else begin

		    legend, [ $
		    string('Night', format='(a-26)')+string(': ',format='(a-5)')+string(night,format='(a-12)'), $
		    string('Prog. ID', format='(a-24)')+string(': ',format='(a-5)')+string(progid,format='(a-13)'), $
		    string('Algorithm', format='(a-24)')+string(': ',format='(a-5)')+string(dumalgo,format='(a-12)'), $
		    string('Field rotation [deg]', format='(a-20)')+string(': ',format='(a-5)')+string(sigfig(fieldrot,4),format='(a-12)'), $
		    string('Obs. time [hr]', format='(a-22)')+string(': ',format='(a-5)')+string(sigfig(tobs,3),format='(a-12)'), $
		    string('DIT science [ms]', format='(a-19)')+string(': ',format='(a-5)')+string(sigfig(dit_sci*1.d3,4),format='(a-12)'), $
		    string('DIT flux [ms]', format='(a-23)')+string(': ',format='(a-5)')+string(sigfig(dit_fl,4),format='(a-12)'), $
		    string('Seeing [asec]', format='(a-21)')+string(': ',format='(a-5)')+string(sigfig(mean(seeing),2),format='(a-12)'), $
		    string(taul+'!D0!N [ms]', format='(a-34)')+string(': ',format='(a-5)')+string(sigfig(mean(tau),3),format='(a-12)'), $
		    string('Wind speed [m/s]', format='(a-18)')+string(': ',format='(a-5)')+string(sigfig(mean(wind),3),format='(a-12)'), $
		    string('RH [%]', format='(a-24)')+string(': ',format='(a-5)')+string(sigfig(mean(rh),3),format='(a-12)'), $
		    string(camfiltdum+' [mag]', format='(a-25)')+string(': ',format='(a-5)')+string(sigfig(mag,4),format='(a-12)')], $
		    box=0, margin=0, /bottom, /left, position=[0.04, 0.09], /normal, charsize=0.7

		  endelse

		device,/close
		set_plot,'x'

		loadct, 0, /silent
		!p.thick=1
		!x.thick=1
		!y.thick=1
		;!p.multi=[0,1,0]

		spawn, '/home/amueller/work/IDLlibs/epstopdf.pl '+fn+'.ps'
        spawn, 'rm '+fn+'.ps'


		;don't even fucking ask....
; 		spawn, 'gs -q -dNOPAUSE -dSAFER -sDEVICE=pdfwrite -sOutputFile='+fn+'.pdf -dMaxSubsetPct=100 -dSubsetFonts=true -dPDFSETTINGS=/prepress -dEmbedAllFonts=true -dAutoRotatePages=/None -r100 -g850x1200 -dPDFFitPage -c "<< /PageOffset [-56 -360] >> setpagedevice" -f '+fn+'.ps -c quit'


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		;---------------------------------------------------------------------------------

  ;if not keyword_set(plotonly) then begin

		;write out detection limits
		if keyword_set(naco) then begin

		  if (agpm eq 1) then begin

		    if (strmatch(datadir, '*Archive*') eq 1) then openw, lun, resdir+'Archive_'+star+'_'+camera_filter[xx]+'_'+yy+mm+dd+'_DetectionLimit'+dirflag+'.txt', width=2000, /get_lun else openw, lun, resdir+star+'_'+camera_filter[xx]+'_'+yy+mm+dd+'_DetectionLimit'+dirflag+'.txt', width=2000, /get_lun

		  endif else begin

		    if (strmatch(datadir, '*Archive*') eq 1) then openw, lun, resdir+'Archive_'+star+'_noAGPM_'+camera_filter[xx]+'_'+yy+mm+dd+'_DetectionLimit'+dirflag+'.txt', /get_lun else openw, lun, resdir+star+'_noAGPM_'+camera_filter[xx]+'_'+yy+mm+dd+'_DetectionLimit'+dirflag+'.txt', width=2000, /get_lun

		  endelse

		endif

		if keyword_set(nici) then begin

		  openw, lun, resdir+star+'_'+camera_filter[xx]+'_'+yy+mm+dd+'_DetectionLimit'+dirflag+'.txt', width=2000, /get_lun

		endif


		if keyword_set(sphere) then openw, lun, resdir+star+'_'+camera_filter[xx]+'_'+yy+mm+dd+'_DetectionLimit'+dirflag+'.txt', width=2000, /get_lun

		idx = where(finite(detection) eq 1)
		contrast = contrast[idx]
		dist = dist[idx]
		detection = detection[idx]
		detection_mag = detection_mag[idx]

		printf, lun, 'Distance_asec DetectionLimit DetectionLimit_MAG'
		for ii=0,n_elements(contrast)-1 do printf, lun, dist[ii], detection[ii], detection_mag[ii], format='(f8.5, f20.15, f11.7)'

		close, lun
		free_lun, lun


		;write out contrast
		if keyword_set(naco) then begin

		  if (agpm eq 1) then begin

		    if (strmatch(datadir, '*Archive*') eq 1) then openw, lun, resdir+'Archive_'+star+'_'+camera_filter[xx]+'_'+yy+mm+dd+'_Contrast'+dirflag+'.txt', width=2000, /get_lun else openw, lun, resdir+star+'_'+camera_filter[xx]+'_'+yy+mm+dd+'_Contrast'+dirflag+'.txt', width=2000, /get_lun

		  endif else begin

		    if (strmatch(datadir, '*Archive*') eq 1) then openw, lun, resdir+'Archive_'+star+'_noAGPM_'+camera_filter[xx]+'_'+yy+mm+dd+'_Contrast'+dirflag+'.txt' else openw, lun, resdir+star+'_noAGPM_'+camera_filter[xx]+'_'+yy+mm+dd+'_Contrast'+dirflag+'.txt', width=2000, /get_lun

		  endelse

		endif

		if keyword_set(nici) then openw, lun, resdir+star+'_'+camera_filter[xx]+'_'+yy+mm+dd+'_Contrast'+dirflag+'.txt', width=2000, /get_lun

		if keyword_set(sphere) then openw, lun, resdir+star+'_'+camera_filter[xx]+'_'+yy+mm+dd+'_Contrast'+dirflag+'.txt', width=2000, /get_lun

		idx = where(finite(detection) eq 1)
		contrast = contrast[idx]
		dist = dist[idx]
		detection = detection[idx]
		detection_mag = detection_mag[idx]

		printf, lun, 'Distance_asec       Contrast Contrast_MAG'
		for ii=0,n_elements(contrast)-1 do printf, lun, dist[ii], contrast[ii], logm(1./contrast[ii]), format='(f8.5, f20.15, f11.7)'

		close, lun
		free_lun, lun

	      endif

; 
; 

	    ;========================================================================

	    endfor

	  endfor

	endfor

      endfor

    endfor

  endfor

endfor

stop
end
