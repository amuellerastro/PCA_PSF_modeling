;;
function rmsmap_am, im, fwhm
  sz = (size(im))[1]
  map = fltarr(sz,sz)
  xc = (sz)/2.0
  yc = (sz)/2.0
  
  map[*]=1.0

  dist_circle, dd, [sz, sz], xc, yc
  
  ;dia=2
  dia = fwhm
  orad=sz/2.0-10
  
  for i=0, orad-1, dia do begin
     w = where(dd gt i and dd le i+dia)
     rms = robust_sigma(im[w])
     ;;rms = stddev(im[w])
     map[w]=rms

  endfor

  return, map
  
end
