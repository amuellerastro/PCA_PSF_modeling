pro fcp_am, fwhm, transm_sci, transm_fl, ditfact, sigma, rad_c, n_br, plsc, dim, camera_filter, right_handed, filtmethod, ADI=adi, PSFFILT=psffilt, $
AGPM=agpm, STUDENT=student, DISPLAY=display, CHECK=check

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

;in case fakes coincide with position of bright CC, offset the branches
offangle = 0.
; offangle = 12.*!dtor


;Read input data
;---------------
obj_tmp = readfits(procdir+'img_'+tg_name_dyn+'_dc.fits',header, /silent)
nref = sxpar(header,'NREF') 

nobj = (size(obj_tmp))[3]; Number of frames in the input object datacube
if (size(obj_tmp))[1] ne dim then MESSAGE, 'Conflicting array size!!!'

;Read the appropriate input PSF
;------------------------------
; r = ROUTINE_INFO('fcp_v83', /SOURCE)
; path = STRMID(r.path,0,STRLEN(r.path)-STRLEN('fcp_v81.pro'))

; if keyword_set(keck)+keyword_set(nici)+keyword_set(naco)+keyword_set(sphere) gt 1 then MESSAGE, 'Please choose only one instrument among Kec, NICI, NACO.'

; if keyword_set('sphere') then begin

  ;IRDIS plate scale is unique and = 12.25 mas/pix

  fcp_orig = readfits(datadir+'PSF_'+tg_name_bas+'.fits', /silent)
  fcp = fcp_orig

;   if (strcompress(camera_filter, /rem) eq 'K1') then fcp=readfits(datadir+'PSF_K1.fits', /silent)
;   if (strcompress(camera_filter, /rem) eq 'K2') then fcp=readfits(datadir+'PSF_K2.fits', /silent)
;   if (strcompress(camera_filter, /rem) eq 'H2') then fcp=readfits(datadir+'PSF_H2.fits', /silent)
;   if (strcompress(camera_filter, /rem) eq 'H3') then fcp=readfits(datadir+'PSF_H3.fits', /silent)
;   if (strcompress(camera_filter, /rem) eq 'J2') then fcp=readfits(datadir+'PSF_J2.fits', /silent)
;   if (strcompress(camera_filter, /rem) eq 'J3') then fcp=readfits(datadir+'PSF_J3.fits', /silent)

  ;in case of IFS YH
;   if (strcompress(camera_filter, /rem) ne 'K1' and strcompress(camera_filter, /rem) ne 'K2' and strcompress(camera_filter, /rem) ne 'H2' and strcompress(camera_filter, /rem) ne 'H3') then begin
;     tmp = readfits(datadir+'PSF.fits',/silent)
;     fcp = tmp[*,*,uint(strcompress(camera_filter, /rem))]
;   endif

; endif

;===========================================================================================

;take DITs and ND filter of science and PSF into account

; fcp = fcp_orig*ditfact*transm_fl/transm_sci

;===========================================================================================

; tmp = fcp
;SO CURRENTLY I AM NOT FILTERING THE PSF!!!!
; if keyword_set(psffilt) then fcp_scaled = mfilter_am(tmp, filtmethod, highpass=2.d0*fwhm, lowpass=0.5d0*fwhm) else fcp_scaled = fcp
; fcp_scaled = fcp
;, /CA)
;===========================================================================================

level_r = readfits(procdir+'vec_'+tg_name_ori+'_level_r.fits', /silent)  ; load the radial noise level file

if keyword_set(check) then begin 

  alpha_r = readfits(procdir+'vec_'+tg_name_bas+'_inter_alpha.fits',header, /silent)
  level_r = level_r/alpha_r

endif

;===========================================================================================

; if keyword_set(display) then begin
; 
;   wset, 2
;   plot, findgen(dim/2)*plsc, sigma*level_r, /ylog, /xlog, yrange=[1e-8,1e-2],xrange=[1e-1,7],xtitle='Angular separation [arcsec]', ytitle='FCP level',linestyle=0, XTicklen=1.0, YTicklen=1.0,XGridStyle=1, YGridStyle=1
; 
; endif

;===========================================================================================


;Compute coordinates and flux level of fake companions
;-----------------------------------------------------
arr_final = dblarr(dim,dim)
fcp_arr = dblarr(dim,dim)
sz_fcp = (size(fcp))[1];/2
if ((double(sz_fcp) mod 2.) eq 0.) then sz_fcp = sz_fcp/2 $
  else sz_fcp = (sz_fcp-1.)/2

;added -1, otherwise mask center shifted by +1,+1 px compared to star center (for dim = even number)
; if ((dim mod 2) ne 0.) then fcp_arr[dim/2-sz_fcp,dim/2-sz_fcp] = fcp $
;     else fcp_arr[dim/2-sz_fcp-1.,dim/2-sz_fcp-1.] = fcp
fcp_arr[dim/2-sz_fcp,dim/2-sz_fcp] = fcp/max(fcp)

n_rad = (size(rad_c))[1] ; Number of different radii where we want to place our fake companions
; coords = intarr(2,n_br,n_rad) ; Coordinates of the center of all the fake companions
;intarr for coords will round the coordinates
coords = dblarr(2,n_br,n_rad) ; Coordinates of the center of all the fake companions

for j = 0, n_br-1 do begin

  theta = j*2*!DPI/n_br + offangle

  for i = 0, n_rad-1 do begin

    r_pix = DOUBLE(rad_c[i])/plsc
    xpos = r_pix * COS(theta)
    ypos = r_pix * SIN(theta)
    coords[0,j,i] = dim/2 + xpos
    coords[1,j,i] = dim/2 + ypos
    ;arr_final = arr_final + (shift_sub(fcp_arr, xpos, ypos) * sigma * level_r[FIX(r_pix)])

    wframe = 150.	;width of frame
    tmp = dblarr(dim+2.*wframe, dim+2.*wframe)
    tmp[wframe:wframe+dim-1, wframe:wframe+dim-1] = fcp_arr
    stmp = fftshift(tmp, xpos, ypos)
    cstmp = stmp[(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1,(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1]
    arr_final = arr_final + (cstmp * sigma * level_r[FIX(r_pix)])
;     arr_final = arr_final + (fftshift(fcp_arr, xpos, ypos) * sigma * level_r[FIX(r_pix)])

  ;window, 0, xs=1000, ys=1000
  ;cgimage, arr_final, /axis, stretch=5
  ;tvcircle, /data, 3., coords[0,j,i], coords[1,j,i], color=cgcolor('yellow')
  
  endfor

  proceeding_text,loop=n_br, i=j, prompt='> Creating FCP Array   '+string(i+1,form='(I4)')

endfor

;Add the fake companions to the input image cube
;---------------------------------------------------
if keyword_set(adi) then begin ; derotate the companion position in case of an ADI sequence

  paral = readfits(datadir+'vec_'+tg_name_bas+'_paral.fits', /silent)

  for i=0,nobj-1 do begin

;     arr_tmp = fltarr(dim,dim)
;     if (right_handed eq 1) then paral[i]=-paral[i]
;     for j = 0, n_br-1 do begin
; 
;       theta = j*2*!DPI/n_br
; 
;       for k = 0, n_rad-1 do begin
; 
;         r_pix = double(rad_c[k])/plsc
;         xpos = r_pix * cos(theta-paral[i]*!DPI/180D0)
;         ypos = r_pix * sin(theta-paral[i]*!DPI/180D0)
; 
; 
; 	wframe = 150.	;width of frame
; 	tmp = dblarr(dim+2.*wframe, dim+2.*wframe)
; 	tmp[wframe:wframe+dim-1, wframe:wframe+dim-1] = fcp_arr
; 	stmp = fftshift(tmp, xpos, ypos)
; 	cstmp = stmp[(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1,(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1]
; 	arr_tmp = arr_tmp + (cstmp * sigma * level_r[FIX(r_pix)])
; 
; ; 	arr_tmp = arr_tmp + (fftshift(fcp_arr, xpos, ypos) * sigma * level_r[FIX(r_pix)])
; 
;       endfor
; 
;     endfor
; 
;     obj_tmp[*,*,i] = obj_tmp[*,*,i] + arr_tmp
    if (right_handed eq 1) then paral[i]=-paral[i]
    obj_tmp[*,*,i] = obj_tmp[*,*,i] + rot(arr_final, paral[i], 1.0, dim/2., dim/2., cubic=-0.5, /pivot)

    proceeding_text,loop=nobj, i=i, prompt='> Add FCP to image cube   '+string(i+1,form='(I4)')

  endfor

endif else begin

  for i=0,nobj-nref-1 do obj_tmp[*,*,i] = obj_tmp[*,*,i] + arr_final

endelse


;Save the new image cube and ancillary data
;------------------------------------------
if keyword_set(check) then begin

  tg_name_dyn = tg_name_dyn+'_fcp_check'
  writefits, procdir+'img_'+tg_name_dyn+'_dc.fits', obj_tmp, header
  writefits, procdir+'img_'+tg_name_bas+'_check_fcparray.fits', arr_final
  writefits, procdir+'vec_'+tg_name_bas+'_check_coords.fits', coords
;   writefits, procdir+'img_'+tg_name_bas+'_check_scaledPSF.fits', fcp_scaled

endif else begin

  tg_name_dyn = tg_name_dyn + '_fcp'
  writefits, procdir+'img_'+tg_name_dyn+'_dc.fits', obj_tmp, header
  writefits, procdir+'img_'+tg_name_bas+'_fcparray.fits', arr_final
  writefits, procdir+'vec_'+tg_name_bas+'_coords.fits', coords
;   writefits, procdir+'img_'+tg_name_bas+'_scaledPSF.fits', fcp_scaled

endelse

end
