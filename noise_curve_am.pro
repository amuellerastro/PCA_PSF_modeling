@showsym.pro

pro noise_curve_am, img, cx, cy, fwhm, plsc, dim, sigma, rin_init, qcut, BOX=box, GAUSSFILT=gaussfilt, $
                        LEVEL_FILE=level_file, ALPHA=alpha, AGPM=agpm, STUDENT=student, DISPLAY=display

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

dim_cc = min([cx,cy,dim-cx,dim-cy]) ; define the minimum FOV radius in the image

img_local = img

; Load self-attenuation file if needed
if keyword_set(alpha) then alpha_r = readfits(procdir+'vec_'+tg_name_bas+'_inter_alpha.fits',header, /silent) else alpha_r = 1d

;========================================================================================

if keyword_set(level_file) then begin ; If the contrast curve has already been computed, use it

  level_r = readfits(procdir+'vec_'+tg_name_ori+'_level_r.fits',header, /silent)  ; load the radial noise level file
  contrast_img = sigma*(level_r)/alpha_r  

endif else begin ; If not, recompute the contrast curve

  ;----------------------------------------------------------------------------------------

  ;cut out planets if needed

  qmode = ''
  if (qcut eq 'y') then begin

    read, 'Cut out or Interpolation of annulus (c/a): ', qmode

    if (dim lt 1000.) then window, 0, xs=1.5*dim, ys=1.5*dim else window, 0, xs=dim, ys=dim
    scale = [min(img_local),max(img_local)]
    cgimage, img_local, /axis, minvalue=scale[0], maxvalue=scale[1]

    quest = ''
    repeat begin

      read, 'Image scale OK? (y/n): ', quest

      if (quest eq 'n') then begin

	print, ''
	print, 'Current values: ', scale
	scale = dblarr(2)
	read, 'Enter new low/high values: ', scale
	cgimage, img_local, minvalue=scale[0], maxvalue=scale[1], /axis

      endif

    endrep until (quest eq 'y')

    !mouse.button = 0
    radius = round(2.*fwhm)
    print, ''
    print, 'Select Candidate(s)'

    i = 0
    xc = dblarr(100) & yc = dblarr(100)

    while !mouse.button ne 4 do begin

      cursor, x, y, 3, /data

      if (!mouse.button eq 1) then begin

	!mouse.button = 0

	x = ceil(x)
	y = ceil(y)

	;refine centering by identifying max intensity
	tmp_cutim = img_local[x-radius:x+radius, y-radius:y+radius]

	mx = max(tmp_cutim, location)
	idxmax = array_indices(tmp_cutim, location)

	xc[i] = x+(idxmax[0]-radius)
	yc[i] = y+(idxmax[1]-radius)

	oplot, xc[i]*[1,1], yc[i]*[1,1], psym=sym(11), color=cgcolor('red'), symsize=2

	i = i+1

      endif

    endwhile

    xsel = xc
    ysel = yc

    idx = where(xsel ne 0.)
    xsel = xsel[idx]
    ysel = ysel[idx]

    for i=0,n_elements(xsel)-1 do begin

      mask_t = shift(dist(dim), xsel[i], ysel[i])
      mask = mask_t ge 2.*fwhm	;CHANGE / ADJUST FACTOR IF NEEDED

      idx = where(mask eq 0)
      idx0 = array_indices(mask,idx)
      ;mask = fix(mask)
      for j=0,n_elements(idx0[0,*])-1 do img_local[idx0[0,j],idx0[1,j]] = !values.f_nan

    endfor

  endif

  ;----------------------------------------------------------------------------------------


;   if keyword_set(gaussfilt) then begin
; 
;     psf = readfits(datadir+'PSF_'+tg_name_bas+'.fits', /silent)
;     psf = psf/max(psf)
;     img_local = convol_fft(img_local, psf)
; 
;   endif

;   if keyword_set(gaussfilt) then img_local = filter_image(img_local, fwhm_gaussian=fwhm/(2.*SQRT(2.*ALOG(2.))))

;   if keyword_set(gaussfilt) then img_local = gauss_smooth(img_local, fwhm/(2.*SQRT(2.*ALOG(2.))), /edge_truncate)

;====================================================================================================

;using apertures at each radii to estimate noise

  r_start = ceil(rin_init*fwhm+(0.5*fwhm))
  r_end = floor(dim_cc-fwhm)
  ndis = dim/2.;r_end-r_start+1
  dis = indgen(ndis);+r_start
  naper = intarr(ndis)	;number of apertures as function of distance from center
  noise = dblarr(ndis)

  for i=r_start,r_end-1 do begin

    naper[i] = (2.*!DPI*dis[i])/fwhm

    ang = 360./naper[i]	;angular sep. between apertures
    theta = dindgen(naper[i])*ang	;angles for each aperture

    xap = dis[i]*cos(theta*!dtor)+dim/2.;+0.5
    yap = dis[i]*sin(theta*!dtor)+dim/2.;+0.5

    ftot = dblarr(naper[i])
    sdev = dblarr(naper[i])
    npix = dblarr(naper[i])

    for j=0,naper[i]-1 do begin

      circint_MOD, img_local, xap[j], yap[j], fwhm, tot, mtot, meantot, maxpx, sdtot, npx, totdif, npxdif, t8

      ftot[j] = meantot
      sdev[j] = sdtot
      npix[j] = npx

    endfor

    ;idx = where(finite(sdev) eq 1)
;     if (qcut eq 'y' and qmode eq 'c') then noise[i] = sqrt((1./(total(npix)))*total(npix*sdev^2.,/nan))
;     if (qcut eq 'y' and qmode eq 'a') then noise[i] = sqrt((1./(total(npix)))*total(npix*sdev^2.))
;     if (qcut eq 'n') then noise[i] = sqrt((1./(total(npix)))*total(npix*sdev^2.))
    noise[i] = stddev(ftot)

;     window, 0, xs=1000, ys=1000
;     cgimage, img_local, /axis, minvalue=min(img_local), maxvalue=max(img_local)
;     for j=0,naper[i]-1 do tvcircle, /data, fwhm/2, xap[j], yap[j], color=cgcolor('green'), thick=2

    proceeding_text, loop=r_end, i=i, prompt='> Noise estimation        '+string(i+1,form='(I4)')

  endfor

  ;fill the noise array with constant values where no full aperture could be set
  noise[0:r_start-1] = noise[r_start]
  noise[r_end:*] = noise[r_end-1]

  if (qcut eq 'y' and qmode eq 'a') then begin

    idxnan = where(finite(noise) ne 1)
    if (idxnan[0] ne -1) then begin

      idx = where(finite(noise) eq 1)
      tmp = interpol(noise[idx], dis[idx], dis)
      noise = tmp

    endif

  endif


  r_fwhm = dindgen(dim_cc)/fwhm
  ; Take into account the effect of the AGPM if needed -- flux loss will be measured so its already taken into account?
;   if keyword_set(agpm) then begin
; 
;     ;equation from D. Mawet, reference Delacroix2013, A&A 553, A98
;     r_fwhm = r_fwhm*0.85 ; correct the effective fwhm
;     coro_thr = 1.d0-exp(-0.831*(r_fwhm)^2) ; transmission profile
; 
;   endif else coro_thr = 1d
  coro_thr = 1d	;correcting this right now with the raw cube

  ; Apply Student's correction for small samples

;   radius = dindgen(dim_cc)
;   penalty = dblarr(dim_cc)
;   nb_resels = fix(round(2.d0*!DPI*radius/fwhm))
;   id_min_radius = (where(radius eq fix(round(fwhm))))[0]
;   false_alarm_p_gaussian =  1d - gauss_pdf(double(sigma))
;   for i=id_min_radius,dim_cc-1 do $
;       penalty[i] = t_cvf(false_alarm_p_gaussian,nb_resels[i]) * sqrt(1.d0+1.d0/(nb_resels[i]-1))/sigma 
  penalty = 1.
  level_r = noise*penalty*(1.d/coro_thr)


;====================================================================================================
;2nd method

;   level_r_th = dblarr(360)
;   level_r = dblarr(dim_cc)
;   tmp = dblarr(360.,2)
; 
;   ;Compute the RMS in boxes if needed
;   if keyword_set(box) then img_local = IMAGE_STDDEV(img_local,box/2D0*fwhm)
; 
;   for rad = 0d, dim_cc-1 do begin
; 
;     ; Transform image in polar coordinates
;     for th = 0d, 360-1 do begin
; 
;       polar_coord = [th,rad]
;       rect_coord = cv_coord(FROM_POLAR=polar_coord, /TO_RECT, /DEGREES) ; convert to polar coordinates to allow azimuthal treatment
;       level_r_th[th] = img_local[round(rect_coord[0])+cx,round(rect_coord[1])+cy]
; 
;     endfor
; 
;     r_fwhm = double(rad)/fwhm
;     ; Take into account the effect of the AGPM if needed
;     if keyword_set(agpm) then begin
; 
;       r_fwhm = r_fwhm*0.85 ; correct the effective fwhm
;       coro_thr = 1.d0-exp(-0.831*(r_fwhm)^2) ; transmission profile
; 
;     endif else coro_thr = 1d
; 
; 
;     ;IF KEYWORD_SET(student) THEN sigma_r = (1.0+(1.3291d0/(r_fwhm-0.64663d0))^1.1595d0) ELSE sigma_r=1d
;     ;sigma_r = (1.0+(1.3291d0/(r_fwhm-0.64663d0))^1.1595d0)	;penalty factor for low number statistics, see Mawet2014
; 
;     ;Compute the noise at the specified radius
; 
;     if keyword_set(box) then level_r[rad] = 1d/coro_thr*median(level_r_th) $
;       else level_r[rad] = 1d/coro_thr*robust_sigma(level_r_th)
;   ;     ELSE level_r[rad] = 1d/coro_thr*sigma_r*ROBUST_SIGMA(level_r_th)
;     ;level_r[rad] = 1d/coro_thr*sigma_r*robust_sigma(level_r_th)
; 
;   endfor
; 
;   ; Apply Student's correction for small samples
; 
;   radius = dindgen(dim_cc)
;   penalty = dblarr(dim_cc)
;   nb_resels = fix(round(2.d0*!DPI*radius/fwhm))
;   id_min_radius = (where(radius eq fix(round(fwhm))))[0]
;   false_alarm_p_gaussian =  1d - gauss_pdf(double(sigma))
;   for i=id_min_radius,dim_cc-1 do $
;       penalty[i] = t_cvf(false_alarm_p_gaussian,nb_resels[i]) * sqrt(1.d0+1.d0/(nb_resels[i]-1))/sigma 
; 
;   level_r = level_r*penalty

;====================================================================================================

  writefits, procdir+'vec_'+tg_name_bas+'_level_r.fits', level_r ; save the radial noise level into a fits file if no noise level file was specified on input
  contrast_img = sigma*(level_r)/alpha_r

endelse


; Plot contrast curve
; -------------------
if keyword_set(display) then begin

  wset, 3
  plot, findgen(dim_cc)*plsc, contrast_img, /ylog, /xlog,xrange=[1e-1,7],xtitle='Angular separation [arcsec]', ytitle='Contrast', charsize=2,$
        linestyle=0, XTicklen=1.0, YTicklen=1.0,XGridStyle=1, YGridStyle=1
  oplot, findgen(dim_cc)*plsc, contrast_img*alpha_r, linestyle=1

endif


; Write fits output
; -----------------

res_fin=fltarr(dim_cc,3)
res_fin[*,0] = findgen(dim_cc)*plsc
res_fin[*,1] = contrast_img
res_fin[*,2] = contrast_img*alpha_r

; if keyword_set(alpha) then begin
; 
;   in = fix(0.15/plsc)
;   out = fix(0.5/plsc)
; 
;   if mean(alpha_r[in:out]) ge 0.1 then begin
; 
; ;     print, '***************************************************************'
; ;     print, '***************************************************************'
; ;     print, 'Target name + reduction params: ', tg_name_dyn 
; ;     print, 'Throughput-corrected X-sigma contrast between 0".15 and 0".5', median(res_fin[in:out,1])
; ;     print, '***************************************************************'
; ;     print, '***************************************************************'
; 
;   endif else begin
; 
;     print, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
;     print, '!!! Signal Attenuation > 10 !!!'
;     print, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
; 
;   endelse
; 
; endif

writefits, procdir+'vec_'+tg_name_dyn+'_contrast.fits', res_fin

END