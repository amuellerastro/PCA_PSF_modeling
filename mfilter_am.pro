;Image filter routine
;D. Mawet, February 2012
;The cutoff inputs must be in pixels

function mfilter_am, image, filtmethod, highpass=highpass, lowpass=lowpass, bandpass_low=bandpass_low, bandpass_high=bandpass_high, median_f=median_f,$ 
mean_f=mean_f, gauss_f=gauss_f, fourier_ideal=fourier_ideal, fourier_gauss=fourier_gauss, fourier_butterworth=fourier_butterworth, ca=ca, display=display

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

if (filtmethod eq 'Wien') then wien = 1
if (filtmethod eq 'Catch') then ca = 1
if (filtmethod eq 'Gaussian Smoothing') then gauss = 1
if (filtmethod eq 'Median Filtering') then median_f = 1
if (filtmethod eq 'Mean Filtering') then mean_f = 1
if (filtmethod eq 'Gaussian Filtering') then gauss_f = 1
if (filtmethod eq 'Fourier') then fourier_ideal = 1
if (filtmethod eq 'Gaussian Fourier') then fourier_gauss = 1
if (filtmethod eq 'Fourier Butterwoth') then fourier_butterworth = 1


print, '*********************************'
print, 'Performing high-pass filtering'
print, '*********************************'

tg_name_dyn=tg_name_dyn+'_filt'

sz = size(image)
if (n_elements(sz) eq 5) then sz3 = 1 else sz3 = sz[3]

for i=0,sz3-1 do begin

  im=image[*,*,i]

  if keyword_set(wien) then begin
  print, 'Doing Wiener filtering for frame #', i+1
      im = im-real_part(fft(fft(im, -1)*wiener(im, /quiet),+1))
      im = filter_image(im, fwhm_gaussian=lowpass)
  endif

  if keyword_set(ca) then begin
  print, 'Doing Catch filtering for frame #', i+1
      im = im-filter_image(im, fwhm_gaussian=highpass)
      im = filter_image(im, fwhm_gaussian=lowpass)
  endif

  if keyword_set(gauss) then begin
  print, 'Doing Gaussian smoothing for frame #', i+1
    im = filter_image(im, fwhm_gaussian=lowpass)
  endif

  if keyword_set(median_f) then begin
  print, 'Doing median filtering for frame #', i+1
	  if keyword_set(highpass) then im=im-filter_image( im, MEDIAN = highpass)
	  if keyword_set(lowpass) then im=filter_image( im, MEDIAN = lowpass)
  endif

  if keyword_set(mean_f) then begin
  print, 'Doing smooth filtering for frame #', i+1
	  if keyword_set(highpass) then im=im-filter_image( im, SMOOTH = highpass)
	  if keyword_set(lowpass) then im=filter_image( im, SMOOTH = lowpass)
  endif

  if keyword_set(gauss_f) then begin
  print, 'Doing Gaussian filtering for frame #', i+1
	  if keyword_set(highpass) then im=im-filter_image( im, FWHM = highpass)
	  if keyword_set(lowpass) then im=filter_image( im, FWHM = lowpass)
  endif


  if keyword_set(fourier_ideal) then begin
  print, 'Doing pure Fourier filtering for frame #', i+1
	  if keyword_set(highpass) then im=bandpass_filter_MOD(im, 1d/highpass, 1d, /ideal)
	  if keyword_set(lowpass) then im=bandpass_filter_MOD(im, 0d, 1d/lowpass, /ideal)
	  if keyword_set(bandpass_low) AND keyword_set(bandpass_high) then $
	  im=bandpass_filter_MOD(im, bandpass_low=bandpass_low, bandpass_high=bandpass_high, /ideal)	
  endif

  if keyword_set(fourier_gauss) then begin
  print, 'Doing Gaussian Fourier filtering for frame #', i+1
	  if keyword_set(highpass) then im=bandpass_filter_MOD(im, 1d/highpass, 1d, /gauss)
	  if keyword_set(lowpass) then im=bandpass_filter_MOD(im, 0d, 1d/lowpass, /gauss)
	  if keyword_set(bandpass_low) AND keyword_set(bandpass_high) then $
	  im=bandpass_filter_MOD(im, bandpass_low=bandpass_low, bandpass_high=bandpass_high, /gauss)	
  endif

  if keyword_set(fourier_butterworth) then begin
  print, 'Doing Butterworth (n=1) Fourier filtering for frame #', i+1
	  if keyword_set(highpass) then im=bandpass_filter_MOD(im, 1d/highpass, 1d, butterworth=1)
	  if keyword_set(lowpass) then im=bandpass_filter_MOD(im, 0d, 1d/lowpass, butterworth=1)
	  if keyword_set(bandpass_low) AND keyword_set(bandpass_high) then $
	  im=bandpass_filter_MOD(im, bandpass_low=bandpass_low, bandpass_high=bandpass_high, butterworth=1)	
  endif

  image[*,*,i]=im

  if keyword_set(display) then begin

    wset, 1
    tvscl, congrid(im,500,500)

  endif

endfor

writefits, procdir+'img_'+tg_name_dyn+'_dc.fits', image, header

return, image

end

