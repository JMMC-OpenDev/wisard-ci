pro wisardgui,input,output,target=target,batch=batch,threshold=threshold,guess=guess,nbiter=nbiter,fov=fov,np_min=np_min,regul=regul,display=display

forward_function WISARD_OIFITS2DATA,WISARD ; GDL does need that, IDL insists on it!

  if (~strmatch(!PATH,'*wisardlib*')) then begin
     b=routine_info('WISARDGUI',/source)
     c=file_dirname(b.path)
     path_to_Wisard=c+"/lib"
     !path=expand_path('+'+path_to_Wisard)+":"+!path
     if (~strmatch(!PATH,'*wisardlib*')) then message,"Unable to set PATH to wisard libraries. Please set up yourself"
  endif


;; define constants
  defsysv, '!DI', dcomplex(0, 1), 1 ;Define i; i^2=-1 as readonly system-variable
  onemas=1d-3*(!DPi/180D)/3600D ; 1 mas in radian

;; passed values
  if (n_elements(threshold) eq 0) then threshold=1d-6                ;  convergence threshold
  if (n_elements(guess) eq 0) then guess=0                       ; default : non guess image
  if (n_elements(nbiter) eq 0) then nbiter=50                     ; number of iterations.
  if (n_elements(fov) eq 0) then fov=14                        ; Field of View of reconstructed image (14*14 marcsec^2 here)
  if (n_elements(np_min) eq 0) then np_min=32                     ;  width of reconstructed image in pixels (same along x and y axis).
  if (n_elements(regul) eq 0) then regul=0 ; gilles code: 0:totvar, 1:psd, 2:l1l2, 3:l1l2_white, 4:support 
  if (n_elements(display) gt 0) then device, decomposed=0, retain=2

  if (n_elements(target) lt 1) then target="*"

  wave_min = -1
  wave_max = -1

;; read input file and get start image & parameters. If no parameters are present (bare oifits) use defaults.
;; FITS_INFO is robust in case of strange tables (no column table). 
fits_info,input, n_ext=next, extname=extname, /silent
extname=strtrim(extname,2)
;read first header
hdu0=mrdfits(input,0,header0)
;if image, take it as guess (however can be updated afterwards, see below)
if n_elements(hdu0) gt 1 then guess=hdu0
;
;examine others
for i=1,next do begin
   if ( extname[i] EQ "IMAGE-OI INPUT PARAM") then begin 
; use a special trick: mrdfits is unable to read this kind of data
; where TFIELDS=0. use fits_read with only headers
      fits_read,input,dummy,inputparamsheader,/header_only,exten_no=i
; parse relevant parameters
      target=strtrim(sxpar(inputparamsheader,"TARGET"),2) ; avoid problems with blanks.
      wave_min=sxpar(inputparamsheader,"WAVE_MIN")
      wave_max=sxpar(inputparamsheader,"WAVE_MAX")
      init_img=strtrim(sxpar(inputparamsheader,"INIT_IMG"),2) ; avoid problems with blanks.
      nbiter=sxpar(inputparamsheader,"MAXITER")
      rgl_name=strtrim(sxpar(inputparamsheader,"RGL_NAME"),2) ; avoid problems with blanks.
      rgl_prio=strtrim(sxpar(inputparamsheader,"RGL_PRIO"),2) ; avoid problems with blanks.
   endif else if ( extname[i] EQ "IMAGE-OI OUTPUT PARAM") then begin ; idem trick
; do nothing
   endif else if extname[i] eq "OI_T3" then begin
      oit3 = ptr_new( mrdfits(input,i,header) )
      oit3head = ptr_new( header )
      oit3arr = (n_elements(oit3arr) gt 0)?[oit3arr,oit3]:oit3
      oit3headarr = (n_elements(oit3headarr) gt 0)?[oit3headarr,oit3head]:oit3head
   endif else if extname[i] eq "OI_VIS2" then begin
      oivis2 = ptr_new( mrdfits(input,i,header) )
      oivis2head = ptr_new( header )
      oivis2arr = (n_elements(oivis2arr) gt 0)?[oivis2arr,oivis2]:oivis2
      oivis2headarr = (n_elements(oivis2headarr) gt 0)?[oivis2headarr,oivis2head]:oivis2head
   endif else if extname[i] eq   "OI_ARRAY" then  begin
      oiarray = ptr_new( mrdfits(input,i,header) )
      oiarrayhead = ptr_new( header )
      oiarrayarr = (n_elements(oiarrayarr) gt 0)?[oiarrayarr,oiarray]:oiarray
      oiarrayheadarr = (n_elements(oiarrayheadarr) gt 0)?[oiarrayheadarr,oiarrayhead]:oiarrayhead
   endif  else if extname[i] eq  "OI_WAVELENGTH" then begin
      oiwave = ptr_new( mrdfits(input,i,header) )
      oiwavehead = ptr_new( header )
      oiwavearr = (n_elements(oiwavearr) gt 0)?[oiwavearr,oiwave]:oiwave
      oiwaveheadarr = (n_elements(oiwaveheadarr) gt 0)?[oiwaveheadarr,oiwavehead]:oiwavehead
   endif else if extname[i] eq  "OI_TARGET" then begin ; only one target.
      oitarget = mrdfits(input,i,targethead)
   endif else begin ; every other tables
      oiother = ptr_new( mrdfits(input,i,header) )
      oiotherhead = ptr_new( header )
      oiotherarr = (n_elements(oiotherarr) gt 0)?[oiotherarr,oiother]:oiother
      oiotherheadarr = (n_elements(oiotherheadarr) gt 0)?[oiotherheadarr,oiotherhead]:oiotherhead
   end
end

; if init_img is defined, then find this hdu, read image and use it as
; guess:
if n_elements(init_img) gt 0 then begin
  w=where(extname eq init_img, count)
  if (count gt 1) then print, 'multiple init images found in input file, discarding.'
  if (count eq 1) then guess=mrdfits(input,w[0],header) ; TODO: use header coordinates etc.
end

;target: if not defined, take first one. If defined, find index, and
;error if target not found:
if (target eq "*") then begin
  target=strtrim(oitarget[0].target,2)
  target_id=oitarget[0].target_id
endif else begin
   targets=strtrim(oitarget.target,2)
   w=where(targets eq target, count)
   if (count lt 1) then message,"target not found in input file."
   target_id=oitarget[w[0]].target_id
endelse

; create table of correspondence: which wave correspond to vis2 and t3
allinsnames=strarr(n_elements(oiwavearr))
for i=0,n_elements(oiwavearr)-1 do allinsnames[i]=sxpar(*oiwaveheadarr[i],"INSNAME")

t3inst=intarr(n_elements(oit3arr))
vis2inst=intarr(n_elements(oivis2arr))

for i=0,n_elements(oivis2arr)-1 do vis2inst[i]=(where(sxpar(*oivis2headarr[i],"INSNAME") eq allinsnames))[0]
for i=0,n_elements(oit3arr)-1 do t3inst[i]=(where(sxpar(*oit3headarr[i],"INSNAME") eq allinsnames))[0]

; all values defined, read data
data=wisard_oifits2data(input,targetname=target)

; simple wavelength selection
if (wave_min lt wave_max) then begin
   w=where(data.wlen ge wave_min and data.wlen le wave_max, count)
   if count lt 1 then message,"no such wavelength range in data"
   data = data[w]
end else begin ; set values of wave_min and max to used values:
   wave_min = min(data.wlen)
   wave_max = max(data.wlen)
end

; start reconstruction. TODO case wrt rgl_name/prior
case regul of
   0: reconstruction = WISARD(data, NBITER = nbiter, threshold=threshold, GUESS = guess, FOV = fov*onemas, NP_MIN = np_min, AUX_OUTPUT = aux_output, /TOTVAR, display=display) 

   1: begin
; if regul == 1 (psd) until PSD keyword is present, create the PSD of
; the documentation:
      distance = double(shift(dist(np_min),np_min/2,np_min/2))
      psd = 1D/((distance^3 > 1D) < 1D6)
      reconstruction = WISARD(data, NBITER = nbiter, threshold=threshold, GUESS = guess, FOV = fov*onemas, NP_MIN = np_min, PSD=psd, AUX_OUTPUT = aux_output, display=display)
   end

   2: begin
      scale = 1D/(NP_min)^2     ; factor for balance between regularization **** and in cost function
      delta = 1D
      reconstruction = WISARD(data, NBITER = nbiter, threshold=threshold, GUESS = guess, FOV = fov*onemas, NP_MIN = np_min, SCALE=scale, DELTA=delta, AUX_OUTPUT = aux_output, display=display)
   end

   3: begin
      scale = 1D/(NP_min)^2     ; factor for balance between regularization **** and in cost function
      delta = 1D
      reconstruction = WISARD(data, NBITER = nbiter, threshold=threshold, GUESS = guess, FOV = fov*onemas, NP_MIN = np_min, SCALE=scale, DELTA=delta, /WHITE, AUX_OUTPUT = aux_output, display=display)
   end
   else: message,"regularization not supported, FIXME!"
endcase
; prepare output HDU
regul_name=['TOTVAR','PSD','L1L2','L1L2WHITE','SOFT_SUPPORT']
FXADDPAR,outhead,'EXTNAME','IMAGE-OI INPUT PARAM'
FXADDPAR,outhead,'TARGET',target
FXADDPAR,outhead,'WAVE_MIN',wave_min
FXADDPAR,outhead,'WAVE_MAX',wave_max
FXADDPAR,outhead,'USE_VIS',0L
FXADDPAR,outhead,'USE_VIS2',1L
FXADDPAR,outhead,'USE_T3',1L
FXADDPAR,outhead,'INIT_IMG','WISARD_IMAGE'
FXADDPAR,outhead,'MAXITER',nbiter
FXADDPAR,outhead,'RGL_NAME',regul_name[regul]
if (n_elements(scale) gt 0 ) then FXADDPAR,outhead,'SCALE',scale
if (n_elements(delta) gt 0 ) then FXADDPAR,outhead,'DELTA',delta
FXADDPAR,outhead,'THRESHOL',threshold
FXADDPAR,outhead,'FOV',fov,'Field of View (mas)'
;
FXADDPAR,imagehead,'EXTNAME','WISARD_IMAGE'
;
; write output file: 
mwrfits,toto,output,header0,/create,/silent,/no_copy,/no_comment
mwrfits,[0],output,outhead,/silent,/no_copy,/no_comment
mwrfits,aux_output.x,output,imagehead,/silent,/no_copy,/no_comment

mwrfits,oitarget,output,targethead,/silent,/no_copy,/no_comment
for i=0,n_elements(oiarrayarr)-1 do mwrfits,*(oiarrayarr[i]),output,*(oiarrayheadarr[i]),/silent,/no_copy,/no_comment
for i=0,n_elements(oiwavearr)-1 do mwrfits,*(oiwavearr[i]),output,*(oiwaveheadarr[i]),/silent,/no_comment ; NO_COPY invalidates use of the table afterwards!
for i=0,n_elements(oiotherarr)-1 do mwrfits,*(oiotherarr[i]),output,*(oiotherheadarr[i]),/silent,/no_copy,/no_comment
; write structures updated with computed values: select only sections
; with target and pass corresponding wavelength.
for i=0,n_elements(oivis2arr)-1 do begin
   col_of_flag=where(strtrim(tag_names(*oivis2arr[i]),2) eq "FLAG", count)
   mwrfits,add_model_oivis2( *oivis2arr[i] , *oiwavearr[vis2inst[i]], aux_output, use_target=target_id ), output, *oivis2headarr[i],/silent,/no_copy,/no_comment,logical_cols=col_of_flag+1
end
for i=0,n_elements(oit3arr)-1 do begin
   col_of_flag=where(strtrim(tag_names(*oit3arr[i]),2) eq "FLAG", count)
   mwrfits,add_model_oit3( *oit3arr[i], *oiwavearr[vis2inst[i]], aux_output, use_target=target_id ), output, *oit3headarr[i],/silent,/no_copy,/no_comment,logical_cols=col_of_flag+1
end
if (n_elements(batch)) then exit
end
