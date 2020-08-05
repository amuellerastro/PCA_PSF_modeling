@readcol.pro
@remchar.pro
@gettok.pro
@strnumber.pro
@min_curve_surf.pro
@cgsurf.pro
@cgsetcolorstate.pro
@cggetcolorstate.pro
@cgcheckforsymbols.pro
@cgdefaultcolor.pro
@colorsareidentical.pro
@cgcolor.pro
@cgsnapshot.pro
@cgcolor24.pro
@cgbitget.pro
@convert_to_type.pro
@t3d.pro
@cgplots.pro
@cgsymcat.pro
@interpol.pro
@sigfig.pro
@cgaxis.pro
@closest.pro
@legend.pro

;currently no correction for SPHERE filters and Johnson filters done

pro planet_mass_detection_limit;, dim, plsc, sigma, hmag, plx, age, filter, sphere=sphere, naco=naco

; common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

;========================================================================

selins = ''
read, 'Select Instrument Naco (n) or Sphere (s): ', selins
; selins = 's'
if (selins eq 'n') then naco = 1
if (selins eq 's') then sphere = 1

;==========================================================================

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
read, 'Select Star?: ', selstar
; selstar = 94
star = stars[selstar-1]
restore, filestar[selstar-1], /verbose

; age = 7.5
; eage = 0.7

if (age le 0.) then begin

  read, 'Enter age in Myr: ', age
  age = age

endif

if (plx le 0. or finite(plx) ne 1) then begin

  read, 'Enter parallax in mas: ', plx
  plx = plx/1.d3

endif

;========================================================================

if keyword_set(sphere) then begin

  datapath = '/home/amueller/work/SPHERE/data/'+star+'/IRDIS/PCA_Reduction/'
  dirs = file_search(datapath+'sPCA*', count=ndir)
  print, ''
  print, 'Data Selection'
  print, '--------------'
  for i=0,ndir-1 do print, strcompress(i+1,/rem), ' ', strmid(dirs[i], strlen(datapath))
  print, ''
  read, 'Select Data: ', sel
; sel = 1
  path = dirs[sel-1]+'/'

  ;------------------------------------------------------------------------

  ;selection of used filter pair
  print, ''
  filter_pair = ['K12', 'H23', 'J23']
  for i=0,n_elements(filter_pair)-1 do print, strcompress(i+1,/rem), ' ', filter_pair[i]
  read, 'Select filter pair: ', selfilt
; selfilt = 2
  if (selfilt eq 1) then filter = ['K1', 'K2']
  if (selfilt eq 2) then filter = ['H2', 'H3']
  if (selfilt eq 3) then filter = ['J2', 'J3']

endif

if keyword_set(naco) then begin

  readcol, '/home/amueller/work/IDLlibs/AO/NACO/datapaths.txt', tmp, format='a', /silent
  print, ''

  idx = where(strmatch(tmp, '*'+star+'*') eq 1)
  if (idx[0] ne -1) then begin

    tmp = tmp[idx]

  endif else begin

    print, ''
    print, 'No data set available. Stop.'
    print, ''
    stop

  endelse


  for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]

  print, ''
  read, 'Select Path: ', selp
;   selp = 1
  datapath = tmp[selp-1]+'PCA_Reduction/'

  if (strmatch(datapath, '*'+star+'*') ne 1) then begin

    print, ''
    print, 'Wrong star selected! Stop.'
    stop

  endif

  dirs = file_search(datapath+'sPCA*', count=ndir)
  print, ''
  print, 'Data Selection'
  print, '--------------'
  for i=0,ndir-1 do print, strcompress(i+1,/rem), ' ', strmid(dirs[i], strlen(datapath))
  print, ''
  read, 'Select Data: ', sel
  path = dirs[sel-1]+'/'

  ;------------------------------------------------------------------------

  filter = ['Lp']

endif

;========================================================================

if keyword_set(sphere) then begin

  if (filter[0] eq 'K1') then begin
    mag = kmag[0]
    emag = ekmag[0]
  endif
  if (filter[0] eq 'J2') then begin
    mag = jmag[0]
    emag = ejmag[0]
  endif
  if (filter[0] eq 'H2') then begin
    mag = hmag[0]
    emag = ehmag[0]
  endif

;   mag = [jmag, hmag, kmag]	;JHK
;   emag = [ejmag, ehmag, ekmag]
endif

if keyword_set(naco) then begin
  mag = lmag[0]	;L'
  emag = elmag[0]	;L'
endif


;==========================================================================



file = file_search(path+'*_Contrast_PCAvar*.txt', count=nf)

if (selins eq 's' and nf eq 2) then begin

  readcol, file[0], dist, t1, t2, format='d,d,d', /silent
  readcol, file[1], t0, t3, t4, format='d,d,d', /silent

  contrast_mag = [[t2],[t4]]

endif

if (selins eq 'n' and nf eq 1) then begin

  readcol, file[0], dist, t1, contrast_mag, format='d,d,d', /silent

endif


;==========================================================================

;check if IR mag is measured by 2MASS

;model
; if keyword_set(sphere) then model_file = '/home/amueller/work/IDLlibs/AO/Models/model.BT-Settl.M-0.0.SPHERE.Vega'
; if keyword_set(naco) then model_file = '/home/amueller/work/IDLlibs/AO/Models/model.BT-Settl.M-0.0.NaCo.Vega'


qmodel = ''
print, ''
print, '1 NextGen (Very Low Mass Stars, no dust formation)'
print, '2 BT-Settl (stars/brown dwarfs/planets (without irradiation), with a cloud model)'
print, '3 AMES-Cond (brown dwarfs/extrasolar planets without irradiation, no dust opacity)'
print, '4 AMES-Dusty (brown dwarfs/extrasolar planets without irradiation, with dust opacity)'
read, 'Select Model: ', qmodel

if keyword_set(sphere) then begin

  if (qmodel eq '1') then model_file = '/home/amueller/work/IDLlibs/AO/Models/model.NextGen.M-0.0.SPHERE.Vega'
  if (qmodel eq '2') then model_file = '/home/amueller/work/IDLlibs/AO/Models/model.BT-Settl.M-0.0.SPHERE.Vega'
  if (qmodel eq '3') then model_file = '/home/amueller/work/IDLlibs/AO/Models/model.AMES-Cond-2000.M-0.0.SPHERE.Vega'
  if (qmodel eq '4') then model_file = '/home/amueller/work/IDLlibs/AO/Models/model.AMES-dusty.M-0.0.SPHERE.Vega'

endif

if keyword_set(naco) then begin

  if (qmodel eq '1') then model_file = '/home/amueller/work/IDLlibs/AO/Models/model.NextGen.M-0.0.NaCo.Vega'
  if (qmodel eq '2') then model_file = '/home/amueller/work/IDLlibs/AO/Models/model.BT-Settl.M-0.0.NaCo.Vega'
  if (qmodel eq '3') then model_file = '/home/amueller/work/IDLlibs/AO/Models/model.AMES-Cond-2000.M-0.0.NaCo.Vega'
  if (qmodel eq '4') then model_file = '/home/amueller/work/IDLlibs/AO/Models/model.AMES-dusty.M-0.0.NaCo.Vega'

endif

if (qmodel eq '1') then mflag = 'NextGen'
if (qmodel eq '2') then mflag = 'BTSettl'
if (qmodel eq '3') then mflag = 'Cond'
if (qmodel eq '4') then mflag = 'Dusty'



;==========================================================================

mj = 1.8986d27
ms = 1.98855d30

;==========================================================================
;==========================================================================
;load model

readcol, model_file, tmp, format='a', delimiter='@', /silent, count=ntot
idx = uintarr(n_elements(tmp))
for i=0,n_elements(tmp)-1 do begin

  flag = strmatch(tmp[i], '*Gyr*')
  if (flag eq 1) then idx[i] = i+1	;+1 to keep the first line

endfor

idx0 = where(idx ne 0)
idx = idx[idx0]-1

nmodels = n_elements(idx)

model_age = dblarr(nmodels)
if keyword_set(sphere) then model = dblarr(9,ntot)	;100 are set as dummy
							;age mass teff J2,J3 / H2,H3 / K1,K2
if keyword_set(naco) then model = dblarr(4,ntot)	;age mass teff Lp

for i=0,nmodels-1 do model_age[i] = double(strmid(tmp[idx[i]], 10, 7))
model_age = model_age*1.d3	;model age in Myr


for i=0,nmodels-1 do begin

  if (i eq 0) then begin

    if keyword_set(sphere) then begin

      nl = idx[i+1]-idx[i]-5
      readcol, model_file, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, format='d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d', numline=nl, skipline=7, /silent

      model[0,0:nl-1] = model_age[i]
      model[1,0:nl-1] = t1
      model[2,0:nl-1] = t2
      model[3,0:nl-1] = t16
      model[4,0:nl-1] = t17
      model[5,0:nl-1] = t13
      model[6,0:nl-1] = t14
      model[7,0:nl-1] = t18
      model[8,0:nl-1] = t19

    endif

    if keyword_set(naco) then begin

      nl = idx[i+1]-idx[i]-5
      readcol, model_file, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, format='d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d', numline=nl, skipline=7, /silent

      model[0,0:nl-1] = model_age[i]
      model[1,0:nl-1] = t1
      model[2,0:nl-1] = t2
      model[3,0:nl-1] = t11

    endif

  endif

  if (i ge 1 and i lt nmodels-2) then begin

    if keyword_set(sphere) then begin

      nl = idx[i+1]-idx[i]-5
      readcol, model_file, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, format='d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d', numline=nl, skipline=idx[i]+7+i*3, /silent

      idx0 = where(model[0,*] eq 0.)
      model[0,idx0[0]:idx0[0]+nl-1] = model_age[i]
      model[1,idx0[0]:idx0[0]+nl-1] = t1
      model[2,idx0[0]:idx0[0]+nl-1] = t2
      model[3,idx0[0]:idx0[0]+nl-1] = t16
      model[4,idx0[0]:idx0[0]+nl-1] = t17
      model[5,idx0[0]:idx0[0]+nl-1] = t13
      model[6,idx0[0]:idx0[0]+nl-1] = t14
      model[7,idx0[0]:idx0[0]+nl-1] = t18
      model[8,idx0[0]:idx0[0]+nl-1] = t19

    endif

    if keyword_set(naco) then begin

      nl = idx[i+1]-idx[i]-5
      readcol, model_file, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, format='d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d', numline=nl, skipline=idx[i]+7+i*3, /silent

      idx0 = where(model[0,*] eq 0.)
      model[0,idx0[0]:idx0[0]+nl-1] = model_age[i]
      model[1,idx0[0]:idx0[0]+nl-1] = t1
      model[2,idx0[0]:idx0[0]+nl-1] = t2
      model[3,idx0[0]:idx0[0]+nl-1] = t11

    endif

  endif

  if (i eq nmodels-1) then begin
 
    if keyword_set(sphere) then begin

      readcol, model_file, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, format='d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d', skipline=idx[i]+7+i*3, count=nl, /silent

      idx0 = where(model[0,*] eq 0.)
      model[0,idx0[0]:idx0[0]+nl-1] = model_age[i]	;age
      model[1,idx0[0]:idx0[0]+nl-1] = t1		;mass
      model[2,idx0[0]:idx0[0]+nl-1] = t2		;Teff
      model[3,idx0[0]:idx0[0]+nl-1] = t16		;J2
      model[4,idx0[0]:idx0[0]+nl-1] = t17		;J3
      model[5,idx0[0]:idx0[0]+nl-1] = t13		;H2
      model[6,idx0[0]:idx0[0]+nl-1] = t14		;H3
      model[7,idx0[0]:idx0[0]+nl-1] = t18		;K1
      model[8,idx0[0]:idx0[0]+nl-1] = t19		;K2

    endif

    if keyword_set(naco) then begin

      readcol, model_file, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, format='d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d', count=nl, skipline=idx[i]+7+i*3, /silent

      idx0 = where(model[0,*] eq 0.)
      model[0,idx0[0]:idx0[0]+nl-1] = model_age[i]	;age
      model[1,idx0[0]:idx0[0]+nl-1] = t1		;mass
      model[2,idx0[0]:idx0[0]+nl-1] = t2		;Teff
      model[3,idx0[0]:idx0[0]+nl-1] = t11		;Lp

    endif

  endif

endfor

idx0 = where(model[0,*] ne 0.)
if (idx0[0] ne -1) then model = model[*,idx0]

model[1,*] = model[1,*]*ms/mj	;convert mass into jupiter masses

;==========================================================================
;==========================================================================

for ff=0,n_elements(filter)-1 do begin

  ;==========================================================================

  ;correction factor between 2MASS and IRDIS NB filters

;   if keyword_set(sphere) then begin

;     mass_filt_path = '/home/amueller/work/IDLlibs/AO/2MASS_RSR/'
;     ird_filt_path = '/home/amueller/work/IDLlibs/AO/transmission_ird/Filters/'
; 
;     if (strmatch(filter[ff], '*J*') eq 1) then begin
; 
;       ;load IRDIS tranmission filter curves
;       readcol, ird_filt_path+'SPHERE_IRDIS_D_J23.txt', li, t1, t2, format='d,d,d', /silent
;       win = [1110., 1345.]
; 
;       ;load 2MASS tranmission filter curves
;       readcol, mass_filt_path+'J_band.txt', lm, fm, format='d,d', /silent
; 
;     endif
; 
;     if (strmatch(filter[ff], '*H*') eq 1) then begin
; 
;       ;load IRDIS tranmission filter curves
;       readcol, ird_filt_path+'SPHERE_IRDIS_D_H23.txt', li, t1, t2, format='d,d,d', /silent
;       win = [1500., 1750.]
; 
;       ;load 2MASS tranmission filter curves
;       readcol, mass_filt_path+'H_band.txt', lm, fm, format='d,d', /silent
; 
;     endif
; 
;     if (strmatch(filter[ff], '*K*') eq 1) then begin
; 
;       ;load IRDIS tranmission filter curves
;       readcol, ird_filt_path+'SPHERE_IRDIS_D_K12.txt', li, t1, t2, format='d,d,d', /silent
;       win = [1990., 2400.]
; 
;       ;load 2MASS tranmission filter curves
;       readcol, mass_filt_path+'K_band.txt', lm, fm, format='d,d', /silent
; 
;     endif
; 
;     if (filter[ff] eq 'J2') then fi = t1
;     if (filter[ff] eq 'J3') then fi = t2
;     if (filter[ff] eq 'H2') then fi = t1
;     if (filter[ff] eq 'H3') then fi = t2
;     if (filter[ff] eq 'K1') then fi = t1
;     if (filter[ff] eq 'K2') then fi = t2
; 
;     lm = lm*1d3
;     lambdainterp = findgen(2500.-900.+1)+900.
;     tm = interpol(fm, lm, lambdainterp)
;     ti = interpol(fi, li, lambdainterp)
;     wl = lambdainterp
; 
;     ;Rayleigh–Jeans law
;     rl = 1.d12/wl^4
; 
;     ;integrate over 2MASS
;     idx = where(tm gt 0.)
;     ftot = total(tm[idx]*rl[idx])
; 
;     ;integrate over both IRDIS filter
;     idx = where(wl gt win[0] and wl lt win[1])
;     ftot_i = total(ti[idx]*rl[idx])
; 
;     ratio = ftot_i/ftot
; 
;     dmag = -2.5d0*alog10(ratio)	;correction factor, i.e. star appears so many magnitudes fainter in IRDIS NB filters
; 
; 
;   endif
; 
;   if keyword_set(naco) then dmag = 0.

  ;==========================================================================


  for xx=0,n_elements(age)-1 do begin

    ;absolute mag of planet
    ;final_mag = contrast_mag[*,ff]+(mag[ff]+dmag)+(5.d0*alog10((plx)+1.))
    final_mag = contrast_mag[*,ff]+(mag)+5.+(5.d0*alog10(plx))	;+dmag

; av ='0'
; final_mag = final_mag-(0.046790767d0*float(av))
;L 0.046790767d0
;H 0.17915115d0
;K 0.11239017d0

    idxl = where(model_age lt age[xx])
    idxr = where(model_age gt age[xx])
    al = model_age[n_elements(idxl)-1]
    ar = model_age[idxr[0]]

    idxc = where(model_age eq age[xx])
    if (idxc[0] ne -1 and eage eq 0.) then begin

      al = age[xx]
      ar = age[xx]

    endif

    if (idxl[0] eq -1) then print, 'Edge reached. Using '+strcompress(round(min(age[xx])),/rem)+' Myr.'
    if (idxr[0] eq -1) then print, 'Edge reached. Using '+strcompress(round(max(age[xx])),/rem)+' Myr.'

    idx = where(model[0,*] ge al and model[0,*] le ar)
    model = model[*,idx]


    if (filter[ff] eq 'J2') then data = reform(model[3,*])
    if (filter[ff] eq 'J3') then data = reform(model[4,*])
    if (filter[ff] eq 'H2') then data = reform(model[5,*])
    if (filter[ff] eq 'H3') then data = reform(model[6,*])
    if (filter[ff] eq 'K1') then data = reform(model[7,*])
    if (filter[ff] eq 'K2') then data = reform(model[8,*])
    if (filter[ff] eq 'Lp') then data = reform(model[3,*])

    if (age[xx] lt 1.) then begin

      print, ''
      print, '******************************************'
      print, 'Age of target younger than youngest model!'
      print, '******************************************'
      print, ''
      stop

    endif

    xout1 = reform(model[0,uniq(model[0,*])])
    tmp = [xout1, age[xx]];, age-eage, age+eage]
    xout = tmp[sort(tmp)]
    xout = xout[uniq(xout)]
    xout = xout[where(xout gt 0.)]

    yout = dindgen(2000.)/2.+0.5	;mass

    mag_int = dblarr(n_elements(xout1),n_elements(yout))
    mass_int = dblarr(n_elements(xout1),n_elements(yout))
    for i=0,n_elements(xout1)-1 do begin

      idx = where(model[0,*] eq xout1[i])
      mag_int[i,*] = interpol(data[idx],model[1,idx],yout)
      mass_int[i,*] = yout

    endfor

    grid = dblarr(n_elements(xout),n_elements(yout))
    if (n_elements(xout) gt 1) then begin

      for i=0,n_elements(yout)-1 do begin

	grid[*,i] = interpol(reform(mag_int[*,i]),xout1,xout)

      endfor

    endif else begin

      grid = mag_int

    endelse

;     xout = dindgen(1000.)/2.+1.d0	;age
;     idx = where(age[xx] le xout)
; 
;     if (idx[0] gt 0) then begin
; 
;       tmp = dblarr(n_elements(xout)+1.)
;       tmp[0:idx[0]-1] = xout[0:idx[0]-1]
;       tmp[idx[0]] = age[xx]
;       tmp[idx[0]+1:*] = xout[idx[0]:*]
; 
;       xout = tmp
; 
;     endif
; 
; 
;     yout = dindgen(2000.)/2.+0.5	;mass
;     x = reform(model[0,*])
;     y = reform(model[1,*])
; 
; ;     x = reform(model[0,*])
; ;     y = reform(model[1,*])
; 
;     if (age gt 300.)  then begin
; 
;       print, ''
;       print, 'Change xout.'
;       stop
; 
;     endif

;     xout = dindgen(300.)+1.d0
;     yout = dindgen(400.)/2.+0.5


;     grid = min_curve_surf(data, x, y, /double, xout=xout, yout=yout)

;     if (xx eq 0) then begin
; 
;       window, 0
;       cgsurf, grid, xout, yout, charsize=4, xtitle='Age [Myr]', ytitle='Mass [M!DJ!N]', ztitle=filter[ff]+' [mag]', /shaded, /save
;       cgplots, x, y, data, psym=2, color='red', /t3d
; 
;     endif

    ;extract model mag of corresponding age

    idx = closest(xout, age[xx])
    age_used = xout[idx[0]]
;     idx = where(xout eq age[xx])
;     if (idx[0] eq -1) then begin
; 
;       print, 'Age does not exist in grid. Stop.'
;       stop
; 
;     endif

    model_mag = reform(grid[idx[0],*])

    ;get mass of corresponding measured mag w.r.t. model mag

    int_mass = interpol(yout, model_mag, final_mag)
    cage = age[xx]

    ;==========================================================================
    ;==========================================================================

  ;   x_interp_out = dindgen(dim/2)
  ;   dist = x_interp_out*plsc


    idx = where(int_mass ge 0.)
    if (idx[0] ne -1) then begin

      dist = dist[idx]
      int_mass = int_mass[idx]

    endif




    if (min(int_mass) lt 0.1) then begin
      print, ''
      print, 'Planet mass below 0.1 M_Jupiter!'
      print, ''
      stop
    endif
    if (min(abs(int_mass)) ge 0.1 and min(abs(int_mass)) lt 1.) then miny = 0.1
    if (min(abs(int_mass)) ge 1. and min(abs(int_mass)) lt 10.) then miny = 1.
    if (min(abs(int_mass)) ge 10. and min(abs(int_mass)) lt 100.) then miny = 10.
    if (min(abs(int_mass)) ge 100. and min(abs(int_mass)) lt 1000.) then miny = 100.

    if (max(int_mass) ge 0.1 and max(int_mass) lt 1.) then maxy = 1.
    if (max(int_mass) ge 1. and max(int_mass) lt 10.) then maxy = 10.
    if (max(int_mass) ge 10. and max(int_mass) lt 100.) then maxy = 100.
    if (max(int_mass) ge 100. and max(int_mass) lt 1000.) then maxy = 1000.

    sigma = 5.

;     set_plot, 'ps'
; 
;     fn = path+'plot_'+filter[ff]+'_contrastMass_'+sigfig(cage,3)+'Myr.ps'
;     device, filename=fn, /color, XSIZE=14, YSIZE=10
;     !p.font=0
;     !p.thick=3
;     !x.thick=3
;     !y.thick=3
; 
;       thisLetter = "163B
;       sigmal = '!9' + String(thisLetter) + '!X'
;       thisLetter = "104B
;       deltal = '!9' + String(thisLetter) + '!X'
; 
;       plot, dist, int_mass, xtitle='Angular separation [arcsec]', ytitle=strcompress(uint(sigma),/rem)+sigmal+' Planet mass detection limit [M!DJ!N]', xr=[0., ceil(max(dist)*10.)/10.], xst=9, yst=1, XTicklen=1.0, YTicklen=1.0, XGridStyle=1, YGridStyle=1, /yl, yr=[miny, maxy], /nodata, position=[0.13,0.11,0.96,0.87]
;       oplot, dist, int_mass, color=cgcolor('blue'), thick=3;, linestyle=2
;       cgAxis, xaxis=1, xrange=[0, ceil(max(dist_au)*10.)/10.], /Save, charsize=1., title='Separation [AU]!C'
;       legend, [star, strcompress(uint(age_used),/rem)+' Myr', sigfig(1./plx,4)+' pc'], box=0, margin=0, /top, /right
; 
;     !p.thick=1
;     !x.thick=1
;     !y.thick=1
; 
;     device,/close
;     set_plot,'x'
;     spawn, 'epstopdf '+fn

    set_plot, 'ps'

    fn = path+'plot_'+filter[ff]+'_contrastMass_'+sigfig(cage,3)+'Myr_'+mflag+'.ps'
    device, filename=fn, /color, XSIZE=14, YSIZE=8.5
    !p.font=0
    !p.thick=3
    !x.thick=3
    !y.thick=3

      thisLetter = "163B
      sigmal = '!9' + String(thisLetter) + '!X'
      thisLetter = "104B
      deltal = '!9' + String(thisLetter) + '!X'

;       cgplot, dist, contrast_mag[*,ff], yr=[ceil(max(contrast_mag[*,ff])),floor(min(contrast_mag[*,ff]))-0.5], yst=9, xtitle='Angular separation [arcsec]', ytitle=strcompress(uint(sigma),/rem)+sigmal+' Contrast ['+deltal+'mag]', xst=9, position=[0.1,0.11,0.87,0.87], charsize=1
; ;       cgaxis, yaxis=1, ylog=1, yr=[1,1000], /save, charsize=1
;       cgAxis, yaxis=1, ylog=1, yr=[1,1000], charsize=1, title=strcompress(uint(sigma),/rem)+sigmal+' Planet mass detection limit [M!DJ!N]', color=cgcolor('blue')
;       cgoplot, dist, int_mass, color=cgcolor('blue')
;       cgAxis, xaxis=1, xrange=[0, ceil(max(dist_au)*10.)/10.], /Save, charsize=1., title='Separation [AU]!C'
;       legend, [star, strcompress(uint(age_used),/rem)+' Myr', sigfig(1./plx,4)+' pc'], box=0, margin=0, /top, /right


      ;dist_au = (1./plx)*206265.d0*tan((dist/3600.)*!dtor)
      dist_au = (1./plx)*206265.d0*tan(([0.1,10.]/3600.)*!dtor)

      yrm = dblarr(2)
      if (min(int_mass) ge 0.1 and min(int_mass) lt 1.) then yrm[0] = 0.1
      if (min(int_mass) ge 1. and min(int_mass) lt 10.) then yrm[0] = 1.
      if (min(int_mass) ge 10. and min(int_mass) lt 100.) then yrm[0] = 10.
      if (min(int_mass) ge 100. and min(int_mass) lt 1000.) then yrm[0] = 100.

      if (max(int_mass) gt 0.1 and max(int_mass) le 1.) then yrm[1] = 1.
      if (max(int_mass) gt 1. and max(int_mass) le 10.) then yrm[1] = 10.
      if (max(int_mass) gt 10. and max(int_mass) le 100.) then yrm[1] = 100.
      if (max(int_mass) gt 100. and max(int_mass) le 1000.) then yrm[1] = 1000.


      cgplot, dist, contrast_mag[*,ff], xr=[0.1,10.], /xl, yr=[ceil(max(contrast_mag[*,ff])),floor(min(contrast_mag[*,ff]))-0.5], yst=9, xtitle='Angular separation [arcsec]', ytitle=strcompress(uint(sigma),/rem)+sigmal+' Contrast ['+deltal+'mag]', xst=9, position=[0.1,0.13,0.87,0.85], charsize=1
;       cgaxis, yaxis=1, ylog=1, yr=[1,1000], /save, charsize=1
      cgAxis, yaxis=1, ylog=1, yr=[yrm[0],yrm[1]], charsize=1, title=strcompress(uint(sigma),/rem)+sigmal+' Planet mass detection limit [M!DJ!N]', color=cgcolor('blue')
      ;cgoplot, dist, int_mass, color=cgcolor('blue')
      cgAxis, xaxis=1, xrange=[min(dist_au),max(dist_au)], /Save, charsize=1., title='Separation [AU]!C',xst=1
      plot, dist, int_mass, /noerase, xst=5, yst=4, position=[0.1,0.11,0.87,0.87], charsize=1, color=cgcolor('blue'), ylog=1, yr=[yrm[0],yrm[1]], xr=[0.1,10.], /xl
      ;legend, [star, strcompress(uint(age_used),/rem)+' Myr', sigfig(1./plx,4)+' pc'], box=0, margin=0, /top, /right

    !p.thick=1
    !x.thick=1
    !y.thick=1

    device,/close
    set_plot,'x'
    spawn, 'epstopdf '+fn


    window, 1
    plot, dist, int_mass, /yl

;     openw, lun, path+star+'_'+filter[ff]+'_contrastMass_'+sigfig(cage,3)+'Myr_'+mflag+'_Av'+av+'.txt', width=1400, /get_lun
    openw, lun, path+star+'_'+filter[ff]+'_contrastMass_'+sigfig(cage,3)+'Myr_'+mflag+'.txt', width=1400, /get_lun

    for i=0,n_elements(dist)-1 do printf, lun, dist[i], int_mass[i], format='(f8.5,f12.4)'


    close, lun
    free_lun, lun

  endfor

endfor

stop
end