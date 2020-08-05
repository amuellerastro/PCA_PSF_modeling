pro candidate_detection, star, camfilt, filt, dim, plsc, fwhm

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

print, 'Candidate detection...'

if (filt eq '1') then resm_fin = mrdfits(procdir+'img_'+tg_name_dyn+'_median.fits', 0, /silent) else resm_fin = mrdfits(procdir+'img_'+tg_name_dyn+'_median_convolved.fits', 0, /silent)

psf = mrdfits(procdir+'PSF_'+camfilt+'.fits',0,/silent)

dimpsf = n_elements(psf[0,*])
if (dimpsf gt 11) then psf = psf[(dimpsf-1)/2.-5:(dimpsf-1)/2.+5,(dimpsf-1)/2.-5:(dimpsf-1)/2.+5]
;psf = psf[20:30,20:30]
psf = psf/max(psf)

image = resm_fin

gauss_noise_std, image, nterms_fit=4, Mode, Std, H, V, Vmean, Hfit
noise = image
noise[*,*] = std

starfinder, image, psf, [5.,5.], 0.8, /rel_threshold, noise_std=noise, n_iter=3, min_distance=fwhm, stars=stars, x, y, fluxes, sigma_x, sigma_y, sigma_f, correlation, /silent

cc = dim/2.
sn_map = resm_fin/rmsmap(resm_fin)

if (total(stars) ne 0.d0) then begin

  ncand = n_elements(x)
  pasf = dblarr(ncand)
  rhosf = dblarr(ncand)
  snsf = dblarr(ncand)

  for i=0,ncand-1 do begin

    rhosf[i] = sqrt((x[i]-cc)^2.+(y[i]-cc)^2.)

    rarad1 = 0.
    rarad2 = (cc-x[i])*plsc
    rarad2 = (rarad2/(3600.*360.))*!dtor
    dcrad1 = 0.
    dcrad2 = (y[i]-cc)*plsc
    dcrad2 = (dcrad2/(3600.*360.))*!dtor
    radif  = rarad2-rarad1
    angle  = atan(sin(radif),cos(dcrad1)*tan(dcrad2)-sin(dcrad1)*cos(radif))
    pasf[i] = angle/!dtor

    if (angle/!dtor lt 0.) then pasf[i] = 360.+angle/!dtor
    if (angle/!dtor gt 360.) then pasf[i] = 360.-angle/!dtor

    snsf[i] = sn_map[x[i],y[i]]

  endfor

  xf = x
  yf = y
  fqsf = correlation

  idx = where(snsf gt 4.5)
  if (idx[0] ne -1) then begin

    ncand = n_elements(idx)

    rhosf = rhosf[idx]
    pasf = pasf[idx]
    xf = xf[idx]
    yf = yf[idx]
    snsf = snsf[idx]
    fqsf = fqsf[idx]

    idxs = bsort(snsf,/reverse)
    rhosf = rhosf[idxs]
    pasf = pasf[idxs]
    xf = xf[idxs]
    yf = yf[idxs]
    snsf = snsf[idxs]
    fqsf = fqsf[idxs]

    rhosf = rhosf*plsc*1.d3


    print, ' #  rho [mas] PA [deg] x [px] y [px]    SNR    Qual.'
    for j=0,ncand-1 do print, j+1, rhosf[j], pasf[j], xf[j], yf[j], snsf[j], fqsf[j], format='(i3,f9.2, f9.1, 2f7.1, f7.1, f9.4)'

  endif else begin

    ncand = 0

    rhosf = !values.d_nan
    pasf = !values.d_nan
    xf = !values.d_nan
    yf = !values.d_nan
    snsf = !values.d_nan
    fqsf = !values.d_nan

  endelse


  fn = procdir+star+'_'+camfilt+'_Candidates.sav'
  save, rhosf, pasf, xf, yf, snsf, fqsf, ncand, filename=fn

endif else begin

  ncand = 0

  rhosf = !values.d_nan
  pasf = !values.d_nan
  xf = !values.d_nan
  yf = !values.d_nan
  snsf = !values.d_nan
  fqsf = !values.d_nan

  fn = procdir+star+'_'+camfilt+'_Candidates.sav'
  save, rhosf, pasf, xf, yf, snsf, fqsf, ncand, filename=fn


endelse


;=========================================================================

; 
; sns_cut = 5.0
; fqs_cut = 0.65
; 
; resm_fin = mrdfits(procdir+'img_'+tg_name_dyn+'_median.fits', 0, /silent)
; 
; ;identify candidates
; 
; sz = dim
; 
; sn_map = resm_fin/rmsmap(resm_fin)
; wc = where(sn_map gt 1.0, numc)
; 
; if (wc[0] ne -1) then begin
; 
;   fqs = fltarr(numc)
;   sns = fltarr(numc)
;   x = fltarr(numc)
;   y = fltarr(numc)
; 
;   for i=0,numc-1 do begin
; 
;     x[i] = wc[i] mod sz
;     y[i] = wc[i] / sz
;     fqs[i] = roundness(sn_map, x[i], y[i], peak=peak, limfw=4, /force)
;     sns[i] = peak
; 
;   endfor
; 
;   ;; reject good points close together
;   wg = where( sns ge sns_cut and fqs le fqs_cut, numg)
;   if (wg[0] ne -1) then begin
; 
;     ch = fltarr(numg)+1
;     test = sn_map
;     xg = x[wg]
;     yg = y[wg]
;     snsg = sns[wg]
;     fqsg = fqs[wg]
;     rhos = fltarr(numg)
;     pas = fltarr(numg)
; 
;     for i=0,numg-1 do begin
; 
;       dist = sqrt((xg[i]-xg)^2 + (yg[i]-yg)^2)
;       ;;wx = where(dist lt 3 and snsg[i] gt snsg)
;       wx = where(dist lt 3)
; 
;       if wx[0] ne -1 then begin
; 
; 	;ch[wx]=0
; 	for j=0,n_elements(wx)-1 do if snsg[i] gt snsg[wx[j]] then ch[wx[j]]=0
; 
;       endif
; 
;     endfor
; 
;     cc = dim/2.
; 
;     ;; best candidates
;     wb = where(ch, numb)
;     test = sn_map
;     for i=0,numb-1 do begin
; 
;       test[ xg[wb[i]],  yg[wb[i]]] += 20
;       rhos[wb[i]] = sqrt( (xg[wb[i]]-cc)^2 + (yg[wb[i]]-cc)^2)
; 
;       rarad1 = 0.
;       rarad2 = (cc-xg[wb[i]])*plsc
;       rarad2 = (rarad2/(3600.*360.))*!dtor
;       dcrad1 = 0.
;       dcrad2 = (yg[wb[i]]-cc)*plsc
;       dcrad2 = (dcrad2/(3600.*360.))*!dtor
;       radif  = rarad2-rarad1
;       angle  = atan(sin(radif),cos(dcrad1)*tan(dcrad2)-sin(dcrad1)*cos(radif))
;       pas[wb[i]] = angle/!dtor
; 
;       if (angle/!dtor lt 0.) then pas[wb[i]] = 360.+angle/!dtor
;       if (angle/!dtor gt 360.) then pas[wb[i]] = 360.-angle/!dtor
; 
;       ;pas[wb[i]] = atan((yg[wb[i]]-cc), (xg[wb[i]]-cc))*180/!pi - 90.0
; 
;     endfor
; 
;     ;; final values
;     xf = xg[wb]
;     yf = yg[wb]
;     snsf = snsg[wb]
;     fqsf = fqsg[wb]
;     rhosf = rhos[wb]
;     pasf = pas[wb]
; 
;     ncand = n_elements(xf)
; 
;     o = bsort(snsf, /reverse)
;     rhosf = rhosf[o]*plsc*1.d3
;     pasf = pasf[o]
;     xf = xf[o]
;     yf = yf[o]
;     snsf = snsf[o]
;     fqsf = fqsf[o]
; 
;     print, ' #  rho [mas] PA [deg] x [px] y [px]    SNR    Qual.'
;     for i=0,ncand-1 do print, i+1, rhosf[i], pasf[i], xf[i], yf[i], snsf[i], fqsf[i], format='(i3,f9.2, f9.1, 2f7.1, f7.1, f9.4)'
; 
;     fn = procdir+star+'_'+camfilt+'_Candidates.sav'
;     save, rhosf, pasf, xf, yf, snsf, fqsf, ncand, sns_cut, fqs_cut, filename=fn
; 
;   endif else begin
; 
;     ncand = 0
; 
;     rhosf = !values.d_nan
;     pasf = !values.d_nan
;     xf = !values.d_nan
;     yf = !values.d_nan
;     snsf = !values.d_nan
;     fqsf = !values.d_nan
; 
;     fn = procdir+star+'_'+camfilt+'_Candidates.sav'
;     save, rhosf, pasf, xf, yf, snsf, fqsf, ncand, sns_cut, fqs_cut, filename=fn
; 
;   endelse
; 
; endif else begin
; 
;   ncand = 0
; 
;   rhosf = !values.d_nan
;   pasf = !values.d_nan
;   xf = !values.d_nan
;   yf = !values.d_nan
;   snsf = !values.d_nan
;   fqsf = !values.d_nan
; 
;   fn = procdir+star+'_'+camfilt+'_Candidates.sav'
;   save, rhosf, pasf, xf, yf, snsf, fqsf, ncand, sns_cut, fqs_cut, filename=fn
; 
; 
; endelse


end