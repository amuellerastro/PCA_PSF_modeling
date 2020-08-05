pro alpha_r_am, fwhm, imgfin_nofcp, n_br, rad_c, dim

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

;test
; spawn, 'cp '+procdir+'img_'+tg_name_bas+'_fcparray.fits '+procdir+'img_'+tg_name_dyn+'_median.fits'
; imgfin_nofcp[*,*] = 0.

imgfin_fcp = READFITS(procdir+'img_'+tg_name_dyn+'_median.fits', /silent)
fcp_array = READFITS(procdir+'img_'+tg_name_bas+'_fcparray.fits', /silent) ;/1000d ; level of

if (size(fcp_array))[1] ne dim then message, 'Conflicting array size!!!' 

fcp_coords = READFITS(procdir+'vec_'+tg_name_bas+'_coords.fits', /silent)
n_rad = (size(fcp_coords))[3]

; Derive the local attenuation
; ----------------------------
alpha_fcp = fltarr(n_br,n_rad) ; Mean value of alpha on the whole fwhm of each fake companion
alpha = fltarr(n_rad) ; Mean azimutal value of alpha for the different radii
alpha_2d = fltarr(dim,dim)

; ratio = median((imgfin_fcp - imgfin_nofcp) / fcp_array, fwhm)
ratio = ((imgfin_fcp - imgfin_nofcp) / fcp_array)

writefits, procdir+'img_'+tg_name_bas+'_ratio.fits', ratio

; Perform "aperture photometry"
; -----------------------------
; anti_aliasing_f = 3
; dim_aa = anti_aliasing_f*dim
; w_t = SHIFT(DIST(dim_aa), dim_aa/2, dim_aa/2)
; w_t = TEMPORARY(w_t) LE CEIL(anti_aliasing_f*1.0*FLOAT(fwhm)/2) ; define an aperture of diameter = FWHM
; w_t = SHIFT(REBIN(DOUBLE(TEMPORARY(w_t)), dim, dim), -dim/2, -dim/2) ; aperture centred on pixel (0,0)

anti_aliasing_f = 1.d0
dim_aa = anti_aliasing_f*dim
; w_t = shift(dist(dim_aa), dim_aa/2-1., dim_aa/2-1.)
w_t = shift(dist(dim_aa), dim_aa/2, dim_aa/2)
w_t = w_t le fwhm
; writefits, procdir+tg_name_bas+'_aperture.fits', w_t

for i=0,n_br-1 do begin

  for j=0,n_rad-1 do begin

    xc = (fcp_coords[*,i,j])[0]
    yc = (fcp_coords[*,i,j])[1]
;     circint_MOD, ratio, xc-1, yc-1, fwhm, tot, mtot, meantot, maxpx, sdtot, npix, totdif, npixdif, t8
    ;alpha_fcp[i,j] = tot/npix

    ;shift the aperture to the companion position
    ;to shift the image using FFT set a frame around the image because of ringing
;here fftshift produces side lobes(?) cannot be used
;     wframe = 150.	;width of frame
;     tmp = dblarr(dim+2.*wframe, dim+2.*wframe)
;     tmp[wframe:wframe+dim-1, wframe:wframe+dim-1] = w_t
;     stmp = fftshift(tmp, xc-dim/2., yc-dim/2.)
;     aper_tmp = stmp[(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1,(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1]

    aper_tmp = shift(w_t, xc-dim/2., yc-dim/2.) ; shift the aperture to the companion position
;     alpha_fcp[i,j] = total((ratio*aper_tmp)[where(aper_tmp gt 0)], /double) / total(aper_tmp,/double)	; aperture "photometry" for the ratio


    tmp = (ratio*aper_tmp);[where(aper_tmp gt 0)]

    idx = where(tmp gt 0.d0)

;     if (idx[0] ne -1) then alpha_fcp[i,j] = robust_mean(tmp[idx]/aper_tmp[idx], 3., sig, numrej) $ ;, goodind=goodind
;       else alpha_fcp[i,j] = robust_mean(tmp/aper_tmp, 3., sig, numrej)	;, goodind=goodind

    ;it happened that there was a large bad pixel(?) which gave alphas>>1
    idxbig = where(tmp[idx] gt 1.)
    if (idxbig[0] ne -1) then begin

      idx1 = where(tmp[idx] le 1.)
      tmpgood = tmp[idx[idx1]]

    endif else begin

      tmpgood = tmp[idx]

    endelse


    if (idx[0] ne -1) then begin

      nboot = 1.d5
      boot_mean, tmpgood, nboot, t1
      alpha_fcp[i,j] = mean(t1)

    endif else begin

      alpha_fcp[i,j] = robust_mean(tmp/aper_tmp, 3., sig, numrej)	;, goodind=goodind

    endelse

; nboot = 1.d4
; boot_mean, tmp[idx], nboot, t1
; at = mean(t1)
; print, at, alpha_fcp[i,j]
; stop


  endfor

  proceeding_text,loop=n_br, i=i, prompt='> Throughput   '+string(i+1,form='(I4)')

endfor

; Azimutal mean of alpha --> alpha(r)
; ------------------------------------
alpha = total(alpha_fcp,1,/double, /nan)/double(n_br)

; Interpolation of alpha(r)
; -------------------------
x_interp_out = indgen(dim/2)
x_interp = reform(fcp_coords[0,0,*]) - dim/2 ; x coordinates of the fake companions along (first) the horizontal branch in the pattern, relative to the image center
interp_alpha = interpol(alpha, x_interp, x_interp_out) ;, /spline)

; Creation of the 2D representation
; ---------------------------------
for x = 0d, dim-1 do begin

  for y = 0d, dim-1 do begin

    rect_coord = [x-dim/2, y-dim/2]
    polar_coord = cv_coord(from_rect=rect_coord, /to_polar, /degrees) 
    if abs(polar_coord[1]) gt dim/2-1 then alpha_2d[x,y] = 0 $
      else alpha_2d[x,y] = interp_alpha[abs(polar_coord[1])]

  endfor

proceeding_text,loop=dim, i=x, prompt='> 2D Map   '+string(i+1,form='(I4)')

endfor


; Final outputs in fits files
; ---------------------------
writefits, procdir+'img_'+tg_name_bas+'_alpha_r_2D.fits', alpha_2d
writefits, procdir+'vec_'+tg_name_bas+'_alpha_r.fits', alpha
writefits, procdir+'vec_'+tg_name_bas+'_inter_alpha.fits', interp_alpha


end