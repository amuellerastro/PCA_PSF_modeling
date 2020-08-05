; PURPOSE:
;   Smart principal component analysis for an ADI image cube
;   Klip and Delta are variable with radius
;   

function SPCA_ADI_VAR, cx_init, cy_init, dim, dim_init, fwhm, plsc, truncate_pca, rin_init, step_init, delta, right_handed, idxrf, flagrf, pointdisk, flagred, candist, VERBOSE=verbose, output=output, section=section

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir


;Rounding and formatting PCA input params
;----------------------------------------
truncate_pca_orig = float(round(truncate_pca))
dim = float(dim)
fwhm_rounded = ceil(fwhm)

; if keyword_set(verbose) then begin
; 
;   print, '***************'
;   print, 'Performing sPCA'
;   print, '***************'
;   print, 'Double-check sPCA params (truncate_pca, delta, rin_init, step_init): ', truncate_pca, delta, rin_init, step_init
; 
; endif


;Reading object input file
;-------------------------
obj = readfits(procdir+'img_'+tg_name_dyn+'_dc.fits',header, /silent)
nobj = (size(obj))[3]
;if (recenter eq 'y') then objrc=readfits(procdir+'img_'+tg_name_ori+'_rsz_shift_dc.fits',headerrc, /silent)

if (strmatch(tg_name_dyn, '*fcp*') ne 1) then im_orig = obj

;Parallactic angle loading
;-------------------------

paral = readfits(datadir+'vec_'+tg_name_bas+'_paral.fits', /silent); reading parallactic angle

;needed for Astrometry_Photometry
if (flagrf eq 1) then paral = paral[idxrf]

; if keyword_set(display) then begin
; 
;   wset, 0
;   plot, paral, xtitle='Frame #', ytitle='Parallactic Angle'
; 
; endif

; Load photometry
; ---------------
; if not keyword_set(normal) then normal = 1
; if normal eq 0 then normal_v = READFITS(procdir+'vec_'+tg_name_bas+'_photometry.fits', /silent) $
;   else normal_v = FLTARR(nobj) + normal  


;Mask center (for coronagraphy or saturated images)
;--------------------------------------------------
;added -1, otherwise mask center shifted by +1,+1 px compared to star center (for dim = even number)
if rin_init ne 0 then begin

;   if ((dim mod 2) ne 0) then mask_t=shift(dist(dim),dim/2,dim/2) $
;     else mask_t=shift(dist(dim),dim/2-1.,dim/2-1.)
  mask_t = shift(dist(dim),dim/2,dim/2)
  mask = mask_t ge fwhm*rin_init
  for i=0,nobj-1 do obj[*,*,i] = obj[*,*,i]*mask

  ;if (recenter eq 'y') then for i=0,nobj-1 do objrc[*,*,i] = objrc[*,*,i]*mask

endif


;Initializing PCA geometrical parameters
;---------------------------------------
objtmp = fltarr(dim*dim,nobj)
mean_obj = fltarr(dim,dim,nobj)
objt = fltarr(dim,dim,nobj)
t = systime(1)
;added -1, otherwise mask center shifted by +1,+1 px compared to star center (for dim = even number)
; if ((dim mod 2) ne 0) then rho = shift(dist(dim),dim/2,dim/2) $; radial coordinate in pixels 
;   else rho = shift(dist(dim),dim/2-1,dim/2-1)
rho = shift(dist(dim),dim/2,dim/2)

rin = rin_init*fwhm ;in pixel, starting at rin_init
rout = dim/2-1 ;finishing at edge of input image
step_law = 1.3
if not keyword_set(section) then n_annuli = floor((rout-rin)/(fwhm*step_init)) $
    else n_annuli = floor((rout-rin)^(1d/step_law)/(step_init*fwhm)) ; total number of radial steps


if keyword_set(section) then begin

  y = replicate(1,dim) # findgen(dim) - (dim/2.)
  y = y / (dim/4)
  y = y*2.*!dpi
  x = transpose(y)
  thetamag = atan(temporary(y),temporary(x))/!dpi*180. + 180.

endif


if (flagred eq 1) then begin

  rin = (candist-(4.*fwhm)) > rin_init*fwhm
  rout = (candist+(5.*fwhm)) < dim/2-1

  n_annuli = floor((rout-rin)/(fwhm*step_init))

endif

;Building zones and performing the PCA
;-------------------------------------
for iann = 0,n_annuli-1 do begin

  if not keyword_set(section) then begin

    in = rin + iann*step_init*fwhm;-dr/2+1 ; inside limit of the zone s in pixel
    step = (rin + (iann+1)*step_init*fwhm) - in ; radial thickness of the zone s in pixel (can be varying)
    out = in + step

  ; 	print, in, out, step
    index_annulus = where( (rho lt out) AND (rho ge in), npix_annulus )
    obj_annulus = fltarr(npix_annulus,nobj) ; create a 2D matrix that will contain the info of the considered annulus for each frame

    if keyword_set(verbose) then print, 'Cutting and saving annulus of '+strcompress(step)+'pix ('+strcompress(step/fwhm)+' fwhm) at a radius of'+strcompress(in)+' pix ('+strcompress(in*plsc)+' arcsec).' 

  endif else begin

    in = rin + floor((iann*step_init*fwhm)^step_law)
    step = (rin + floor(((iann+1)*step_init*fwhm)^step_law)) - in
    out = in+step

    g = 0.75; geometry of the optimization zone (see Lafreniere et al. 2007, ApJ 660)
		  ; 0.75* = elongated along theta, squashed radially
		  ; 1     = sqr
		  ; 2     = elongated radially, 2xlong
    Na = 500.;250.
    fc = 360.
    dphi = (g/2+(2*in/fwhm)*sqrt(g/(!dpi*Na)))^(-1d)
    dphi = (dphi/!dpi*180)

    nq = floor(fc/dphi)
    if (nq mod 2) eq 0 then nq = nq+1

    if keyword_set(verbose) then print, 'Cutting and saving annulus of '+strcompress(step)+'pix ('+strcompress(step/fwhm)+' fwhm) at a radius of'+strcompress(in)+' pix ('+strcompress(in*plsc)+' arcsec).' 

  endelse


  if not keyword_set(section) then begin

    ;Fill the 2D matrix with the pixels intensities from the considered annulus
    for fr=0, nobj-1 do obj_annulus[*,fr] = (obj[*,*,fr])[index_annulus]

    ;Conditioning data (removing mean, necessary for PCA)
    ;if keyword_set(verbose) then print, 'Doing sPCA on annulus '+strcompress(iann+1)+' of '+strcompress(n_annuli)
    data=0 & covMatrix=0 & eigenval=0 & eigenvect=0 
    obj_annulus_mean=fltarr(nobj)
  
    for fr=0,nobj-1 do begin

  ;     obj_annulus_mean[fr] = mean(obj_annulus[*,fr])
      obj_annulus_mean[fr] = median(obj_annulus[*,fr])
      obj_annulus[*,fr] = obj_annulus[*,fr]-obj_annulus_mean[fr]

    endfor


    ;Doing the sPCA analysis frame by frame
    for fr=0,nobj-1 do begin 

      ;===================================================================

      ;maximise exclusion angle where we still have enough frames and are able to use e.g. klip=5
      if (pointdisk eq 'd') then kliprange = truncate_pca_orig else kliprange = dindgen(truncate_pca_orig-4.)+5.

      ; delta_v = linspace(0.1,1.5,29)
      delta_v = linspace(0.01,delta,19)
      ;delta_v = linspace(0.1,0.75,14)
      delta_v_deg = delta_v / (in/fwhm) / !dpi*180.


      nref_v = lonarr(n_elements(delta_v))
      for i=0,n_elements(delta_v)-1 do nref_v[i] = n_elements(where((paral lt paral[fr]-delta_v_deg[i]) or (paral gt paral[fr]+delta_v_deg[i])))


      nref_v2 = lonarr(n_elements(nref_v),n_elements(kliprange))
      for i=0,n_elements(kliprange)-1 do nref_v2[*,i] = nref_v-kliprange[i]

      ;find the klip value for the largest possible angle
      idx = where(nref_v2 ge 0)
      if (idx[0] eq -1) then begin

	print, ''
	print, 'No unique parameter set found. Stop.'
	print
	stop

      endif


      i = n_elements(nref_v)
      j = n_elements(kliprange)

      if (pointdisk eq 'p') then begin

	for xx=i-1,0,-1 do begin

	  for yy=j-1,0,-1 do begin

	    ;print, xx, yy, nref_v2[xx,yy]
	    if (xx lt 0 or yy lt 0) then stop

	    if (nref_v2[xx,yy] ge 0) then goto, jump1

	  endfor

	endfor

  jump1:

	i = xx
	j = yy

      endif else begin

	repeat begin

	  i = i-1

	endrep until (nref_v2[i] ge 0)

      endelse


      delta_v = delta_v[i]
      if (pointdisk eq 'p') then truncate_pca = kliprange[j]
      delta_deg = delta_v_deg[i]

      if keyword_set(output) then begin

	openw, lun, procdir+'Annulus_'+strcompress(iann+1,/rem)+'_KLIP_Delta.txt', width=1400, /get_lun, /append
	  ;if (iann eq 0 and fr eq 0) then printf, lun, 'Frame Delta KLIP'
	  printf, lun, fr+1, delta_v, truncate_pca, format='(i5, f6.2, i5)'
	close, lun
	free_lun, lun

      endif

      ;===================================================================

      index_ref = where((paral lt paral[fr]-delta_deg) or (paral gt paral[fr]+delta_deg)) ; pick out the frames that are not in the parallactic exclusion zone 
      nref = n_elements(index_ref)
      ;print, strcompress(iann+1)+'/'+strcompress(n_annuli), fr, nref
      data = obj_annulus[*,index_ref]

; print, truncate_pca, paral[fr]-delta_deg, paral[fr]+delta_deg
; hak
      ;PCA	
      covMatrix = matrix_multiply(data, data, /ATRANSPOSE)

      eigenval = LA_EIGENQL(covMatrix, EIGENVECTORS=eigenvect, range=[nref-truncate_pca,nref-1], /double)
      eigenval = reverse(eigenval)
      eigenvect = reverse(eigenvect,2)


      ;if keyword_set(verbose) then print, 'PCA Eigenvalues: ', eigenval

      ;Principal components 
      pc = matrix_multiply(eigenvect, data, /ATRANSPOSE, /BTRANSPOSE)
      supref = fltarr(npix_annulus)

      for k=0,truncate_pca-1 do begin

	pc[k,*] = pc[k,*] / sqrt(eigenval[k])
	coeff = total(obj_annulus[*,fr]*reform(pc[k,*]))
	;print, coeff
	supref = supref+coeff*reform(pc[k,*])

      endfor

  ;     if (recenter eq 'y') then objtmp[index_annulus,fr] = obj_annulus_obj[*,fr]-supref $
  ;       else objtmp[index_annulus,fr] = obj_annulus[*,fr]-supref;+obj_annulus_mean[fr]
      objtmp[index_annulus,fr] = obj_annulus[*,fr]-supref;+obj_annulus_mean[fr]

    endfor	;nobj

  endif	;annulus only

  if keyword_set(section) then begin

    for q=0,nq-1 do begin	;going through section of each annulus

      index_annulus = where((thetamag ge (q*fc/nq)) and (thetamag le ((q+1)*fc/nq)) and (rho lt out) and (rho ge in), npix_annulus)
      obj_annulus = fltarr(npix_annulus,nobj) ; create a 2D matrix that will contain the info of the considered 	section of the annulus for each frame

      for fr=0,nobj-1 do obj_annulus[*,fr] = (obj[*,*,fr])[index_annulus]

      data = 0 & covMatrix = 0 & eigenval = 0 & eigenvect = 0 
      obj_annulus_mean = fltarr(nobj)
    
      for fr=0,nobj-1 do begin

        ;obj_annulus_mean[fr] = mean(obj_annulus[*,fr])
	obj_annulus_mean[fr] = median(obj_annulus[*,fr])
	obj_annulus[*,fr] = obj_annulus[*,fr]-obj_annulus_mean[fr]

      endfor

      ;Doing the sPCA analysis frame by frame
      for fr=0,nobj-1 do begin 

	;===================================================================

	;maximise exclusion angle where we still have enough frames and are able to use e.g. klip=5
	if (pointdisk eq 'd') then kliprange = truncate_pca_orig else kliprange = dindgen(truncate_pca_orig-4.)+5.
	; delta_v = linspace(0.1,1.5,29)
	delta_v = linspace(0.01,delta,19)
	;delta_v = linspace(0.1,0.75,14)
	delta_v_deg = delta_v / (in/fwhm) / !dpi*180.

	nref_v = lonarr(n_elements(delta_v))
	for i=0,n_elements(delta_v)-1 do nref_v[i] = n_elements(where((paral lt paral[fr]-delta_v_deg[i]) or (paral gt paral[fr]+delta_v_deg[i])))


	nref_v2 = lonarr(n_elements(nref_v),n_elements(kliprange))
	for i=0,n_elements(kliprange)-1 do nref_v2[*,i] = nref_v-kliprange[i]

	;find the klip value for the largest possible angle
	idx = where(nref_v2 ge 0)
	if (idx[0] eq -1) then begin

	  print, ''
	  print, 'No unique parameter set found. Stop.'
	  print
	  stop

	endif


	i = n_elements(nref_v)
	j = n_elements(kliprange)

	if (pointdisk eq 'p') then begin

	  for xx=i-1,0,-1 do begin

	    for yy=j-1,0,-1 do begin

	      ;print, xx, yy, nref_v2[xx,yy]
	      if (xx lt 0 or yy lt 0) then stop

	      if (nref_v2[xx,yy] ge 0) then goto, jump2

	    endfor

	  endfor

    jump2:

	  i = xx
	  j = yy

	endif else begin

	  repeat begin

	    i = i-1

	  endrep until (nref_v2[i] ge 0)

	endelse


	delta_v = delta_v[i]
	if (pointdisk eq 'p') then truncate_pca = kliprange[j]
	delta_deg = delta_v_deg[i]

	if keyword_set(output) then begin

	  openw, lun, procdir+'Annulus_'+strcompress(iann+1,/rem)+'_KLIP_Delta.txt', width=1400, /get_lun, /append
	    ;if (iann eq 0 and fr eq 0) then printf, lun, 'Frame Delta KLIP'
	    printf, lun, fr+1, delta_v, truncate_pca, format='(i5, f6.2, i5)'
	  close, lun
	  free_lun, lun

	endif

	;===================================================================


	index_ref = where((paral lt paral[fr]-delta_deg) or (paral gt paral[fr]+delta_deg)) ; pick out the frames that are not in the parallactic exclusion zone 
	nref = n_elements(index_ref)
	;print, strcompress(iann+1)+'/'+strcompress(n_annuli), fr, nref
	data = obj_annulus[*,index_ref]

	;PCA	
	covMatrix = matrix_multiply(data, data, /ATRANSPOSE)
	eigenval = LA_EIGENQL(covMatrix, EIGENVECTORS=eigenvect, range=[nref-truncate_pca,nref-1], /double)
	eigenval = reverse(eigenval)
	eigenvect = reverse(eigenvect,2)

	;if keyword_set(verbose) then print, 'PCA Eigenvalues: ', eigenval

	;Principal components 
	pc = matrix_multiply(eigenvect, data, /ATRANSPOSE, /BTRANSPOSE)
	supref = fltarr(npix_annulus)	

	for k=0,truncate_pca-1 do begin

	  pc[k,*] = pc[k,*] / sqrt(eigenval[k])
	  coeff = total(obj_annulus[*,fr]*reform(pc[k,*]))
	  ;print, coeff
	  supref = supref+coeff*reform(pc[k,*])

	endfor

    ;     if (recenter eq 'y') then objtmp[index_annulus,fr] = obj_annulus_obj[*,fr]-supref $
    ;       else objtmp[index_annulus,fr] = obj_annulus[*,fr]-supref;+obj_annulus_mean[fr]
	objtmp[index_annulus,fr] = obj_annulus[*,fr]-supref;+obj_annulus_mean[fr]

      endfor	;nobj

    endfor

  endif	;section in annulus

endfor	;annulus

if keyword_set(verbose) then print, 'Timing for annulus sPCA: ', systime(1)-t

objt = reform(objtmp,dim,dim,nobj)

objt_norot = objt
var = variance(objt_norot, dim=3, /double, /nan)
var_cube = cmreplicate(var,n_elements(objt[0,0,*]))


if (strmatch(tg_name_dyn, '*fcp*') ne 1) then imrot = im_orig

;Reconstructing and derotating
;-----------------------------
if keyword_set(verbose) then print, 'Field rotation [deg]: ', abs(paral[0]-paral[n_elements(paral)-1])

for i=0, nobj-1 do begin

  ;if keyword_set(verbose) then print, 'Derotate by: ', paral[i]
  if right_handed eq 1 then paral[i]=-paral[i]

; ;added -1, otherwise mask center shifted by +1,+1 px compared to star center (for dim = even number)
;   if ((dim mod 2) ne 0) then objt[*,*,i]=rot(objt[*,*,i],-paral[i],1.0,dim/2,dim/2,cubic=-0.5,/pivot) $
;     else objt[*,*,i]=rot(objt[*,*,i],-paral[i],1.0,dim/2-1.,dim/2-1.,cubic=-0.5,/pivot)
;added -1, otherwise mask center shifted by +1,+1 px compared to star center (for dim = even number)
  objt[*,*,i] = rot(objt[*,*,i],-paral[i],1.0,dim/2,dim/2,cubic=-0.5,/pivot)

  var_cube[*,*,i] = rot(var_cube[*,*,i],-paral[i],1.0,dim/2,dim/2,cubic=-0.5,/pivot)

  if (strmatch(tg_name_dyn, '*fcp*') ne 1) then imrot[*,*,i]=rot(im_orig[*,*,i],-paral[i],1.0,dim/2,dim/2,cubic=-0.5,/pivot)

endfor

out_name = tg_name_dyn
tg_name_dyn=tg_name_dyn+'_pca'


;Final images and outputs
;------------------------

resm_fin = median(objt,dim=3)

; resm_fin= filter_image(resm_fin, fwhm_gaussian=0.5*fwhm)
writefits, procdir+'img_'+tg_name_dyn+'_median.fits',resm_fin, header

;Bottom 2017, not using it, many artifacts in final image, no much of an improvement so far compared to median
i1 = 1./(total(1./var_cube, 3, /double, /nan))
i2 = total(objt/var_cube,3, /double, /nan)
final_image = i1*i2
writefits, procdir+'img_'+tg_name_dyn+'_Bottom17.fits', final_image, header

sd = stddev(objt, dimension=3)
;sem = sd ;1.253d0*sd/sqrt(nobj)
sem = dblarr(dim,dim)
for i=0,dim-1 do begin

  for j=0,dim-1 do begin

    sem[i,j] = median(abs(objt[i,j,*]-median(objt[i,j,*])))

  endfor

endfor

writefits, procdir+'img_'+tg_name_dyn+'_median_sem.fits',sem, header

writefits, procdir+tg_name_dyn+'_dc.fits', objt

if (strmatch(tg_name_dyn, '*fcp*') ne 1) then begin

  writefits, procdir+'img_'+out_name+'_median.fits',  median(imrot,dim=3), header
  sd = stddev(imrot, dimension=3)
  sem = sd ;1.253d0*sd/sqrt(nobj)
  writefits, procdir+'img_'+out_name+'_median_sem.fits',sem, header

endif


return, resm_fin

end
