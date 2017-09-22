;
;+
;
; NAME: WISARDGUI
;
; PURPOSE: 
;     fronted of wisard interferometric image reconstruction
;     procedure.
;       
; CATEGORY:
;     Inverse Problems
;
; CALLING SEQUENCE: 
;     WISARDGUI,INPUT,OUTPUT,
;     TARGET=TARGET,INTERACTIVE=INTERACTIVE,THRESHOLD=THRESHOLD,
;     GUESS=GUESS,NBITER=NBITER,FOV=FOV,NP_MIN=NP_MIN,REGUL=REGUL,
;     POSITIVITY=POSITIVITY,OVERSAMPLING=OVERSAMPLING,
;     INIT_IMG=INIT_IMG,RGL_PRIO=RGL_PRIO,DISPLAY=DISPLAY,
;     MU_SUPPORT=MU_SUPPORT, FWHM=FWHM, WAVERANGE=WAVERANGE,
;     SIMULATED_DATA=SIMULATED_DATA, USE_FLAGGED_DATA=USE_FLAGGED_DATA
;
; KEYWORD PARAMETERS:
;    INPUT: input OIFITS or (better) OImaging OIFITS file (contains
;    already guess image and parameters. Otherwise, use Keywords to enter
;    these parameters
;    OUTPUT: product  OImaging OIFITS file name.
;
;    TARGET: the object name (in case there are many in the input
;    file)
;    NBITER: max number of iterations (50 by default)
;    NP_MIN: minimum number of resels to reconstruct. default computed
;    internally 
;    GUESS: Guess image (fits, seed doc for header parameters). Not mandatory.
;    INIT_IMG: start with this image. Very useful for a start.
;
;    The other kw are best described in the documentation, see:
;    WISARD:
;    http://www.mariotti.fr/doc/approved/JMMC-MAN-2500-0001.pdf
;    OImaging interface:
;    http://www.mariotti.fr/doc/approved/JRA4Fp7Report2016.pdf
;
; LICENCE:
; Copyright (C) 2017, G. Duvert unless explicitly mentioned in files.
;
; This program is free software; you can redistribute it and/or modify  
; it under the terms of the GNU General Public License as published by  
; the Free Software Foundation; either version 2 of the License, or     
; (at your option) any later version.                                   
;
;-

pro wisardgui,input,output,target=target,interactive=interactive,threshold=threshold,guess=guess,nbiter=nbiter,fov=fov,np_min=np_min,regul=regul,positivity=positivity,oversampling=oversampling,init_img=init_img,rgl_prio=rgl_prio,display=display,mu_support=mu_support, fwhm=fwhm, waverange=waverange, simulated_data=issim, use_flagged_data=use_flagged_data

  forward_function WISARD_OIFITS2DATA,WISARD ; GDL does not need that, IDL insists on it!

  wisard_ci_version=getenv('WISARD_CI_VERSION')
  if (strlen( wisard_ci_version) lt 1) then wisard_ci_version='Unversioned'
  message,/informational,"Welcome to Wisard, version "+wisard_ci_version+", you have accepted the copyrights."

  dotarget=n_elements(target) ne 0
  dothreshold=n_elements(threshold) ne 0
  doguess=n_elements(guess) ne 0
  donbiter=n_elements(nbiter) ne 0
  dofov=n_elements(fov) ne 0
  donp_min=n_elements(np_min) ne 0
  doregul=n_elements(regul) ne 0
  doovers=n_elements(oversampling) ne 0
  doposit=n_elements(positivity) ne 0
  doinit_img=n_elements(init_img) ne 0
  dorgl_prio=n_elements(rgl_prio) ne 0
  dowaverange=n_elements(waverange) ne 0
  
;; define constants
  regul_name=['TOTVAR','PSD','L1L2','L1L2WHITE','SOFT_SUPPORT']
  defsysv, '!DI', dcomplex(0, 1), 1 ;Define i; i^2=-1 as readonly system-variable
  onemas=1d-3*(!DPi/180D)/3600D     ; 1 mas in radian
  
;; defaults in absence of passed values. Normally these should be read in one of the headers.
  if ~dotarget    then target="*"
  if ~dothreshold then threshold=1d-6 else threshold=double(threshold); convergence threshold
  if ~doguess     then guess=0        ; default : non guess image
  if ~donbiter    then nbiter=50      else nbiter=fix(nbiter) ; number of iterations.
  if ~dofov       then fov=14.0D      else fov = double(fov) ; Field of View of reconstructed image (14*14 marcsec^2 here)
  if ~donp_min    then np_min=32L      else np_min = fix(np_min) ; width of reconstructed image in pixels (same along x and y axis).
  if ~doregul     then regul=0  else begin ; TOTVAR by default
      w=where(strtrim(regul_name,2) eq regul, count) ; regul is a name here
      if count lt 1 then regul=0 else regul=w[0]       ; regul is now an integer: 0:totvar, 1:psd, 2:l1l2, 3:l1l2_white, 4:support 
  endelse
  if ~doovers     then oversampling=1 else oversampling = fix (oversampling) ; oversampling
  if ~doposit     then positivity=1   else positivity = fix (positivity) ; quoted: "misplaced curiosity"

  if n_elements(display) gt 0 then device, decomposed=0, retain=2
  if n_elements(display) gt 0 then interactive=1 

; catch any error to exit cleanly in batch mode (GDL does not - yet)
  if ~n_elements(interactive) then begin
     CATCH, error_status
     if (error_status ne 0) then begin
        CATCH,/CANCEL
        message,/reissue_last,/informational
        print,"error has occured, exiting."
        exit, status=1
     end
  end

; after catch to handle this first error message if necessary. 
  if (~strmatch(!PATH,'*wisardlib*')) then begin
     b=routine_info('WISARDGUI',/source)
     c=file_dirname(b.path)
     path_to_Wisard=c+"/lib"
     !path=expand_path('+'+path_to_Wisard)+":"+!path
     if (~strmatch(!PATH,'*wisardlib*')) then message,"Unable to set PATH to wisard libraries. Please set up yourself"
  endif

  wave_min = -1
  wave_max = -1
  if (dowaverange) then begin
     wave_min=waverange[0]*1D-6
     wave_max=waverange[1]*1D-6
  endif
;; read input file and get start image & parameters. If no parameters are present (bare oifits) use defaults.
;; FITS_INFO is robust in case of strange tables (no column table). 
  fits_info,input, n_ext=next, extname=extname, /silent
  if (next le 4) then message,"Input file is probably not an OIFITS file at all!"
  extname=strtrim(extname,2)
;;create list of hdu numbers; and list of hdunames(if any) for further use:
  hdunames=replicate('',next)
  hdunumbers=indgen(next)

;read first header
  hdu0=mrdfits(input,0,main_header)
;if image, take it as guess (however can be updated afterwards, see below)
  if n_elements(hdu0) gt 1 then guess=hdu0
; eventually, if has an HDUNAME, get it
  hduname=strtrim(sxpar(main_header,"HDUNAME"),2)
  if hduname ne '0' then hdunames[0]=hduname
;
;examine others
  for i=1,next do begin
     if ( extname[i] EQ "IMAGE-OI INPUT PARAM") then begin 
; use a special trick: mrdfits is unable to read this kind of data
; where TFIELDS=0. use fits_read with only headers
        fits_read,input,dummy,inputparamsheader,/header_only,exten_no=i
; parse relevant parameters
        if ~dotarget then begin
           req_target=strtrim(sxpar(inputparamsheader,"TARGET"),2) ; avoid problems with blanks.
           if req_target ne '0' then target=req_target
        end
        if ~donbiter then begin
           req_nbiter=sxpar(inputparamsheader,"MAXITER")
           if req_nbiter gt 0 then nbiter=req_nbiter
        end
        if ~doregul then begin
           req_rgl_name=strtrim(sxpar(inputparamsheader,"RGL_NAME"),2) ; avoid problems with blanks.
           if req_rgl_name ne '0' then begin
              w=where(strtrim(req_rgl_name,2) eq regul_name, count)
              if count lt 1 then regul=0 else regul=w[0] ; avoid : ;  message,"unrecognized regularisation name."
           end 
        end
        if ~dofov then begin
           req_fov=sxpar(inputparamsheader,"FOV") 
           if req_fov gt 0 then fov=req_fov
        end
        if ~dothreshold then begin
           req_threshold=sxpar(inputparamsheader,"THRESHOL") 
           if req_threshold gt 0 then threshold=req_threshold
        end
        if ~donp_min then begin
           req_np_min=sxpar(inputparamsheader,"NP_MIN") 
           if req_np_min gt 0 then np_min=req_np_min
        end
        if ~doovers then begin
           req_overs=sxpar(inputparamsheader,"OVERSAMP") 
           if req_overs gt 0 then oversampling=req_overs
        end
        if ~doinit_img then init_img=strtrim(sxpar(inputparamsheader,"INIT_IMG"),2) ; avoid problems with blanks.
        if ~dorgl_prio then rgl_prio=strtrim(sxpar(inputparamsheader,"RGL_PRIO"),2) ; avoid problems with blanks.

        if ~dowaverange then wave_min=sxpar(inputparamsheader,"WAVE_MIN")
        if ~dowaverange then wave_max=sxpar(inputparamsheader,"WAVE_MAX")
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
     endif else begin           ; every other tables: may contain an hduname image
        oiother = ptr_new( mrdfits(input,i,header) )
        ; eventually, if has an HDUNAME, get it
        hduname=strtrim(sxpar(header,"HDUNAME"),2)
        if hduname ne '0' then hdunames[i]=hduname
        oiotherhead = ptr_new( header )
        oiotherarr = (n_elements(oiotherarr) gt 0)?[oiotherarr,oiother]:oiother
        oiotherheadarr = (n_elements(oiotherheadarr) gt 0)?[oiotherheadarr,oiotherhead]:oiotherhead
     end
  end

; if init_img is defined in header, then find this hdu, read image and use it as guess
  if ~doinit_img then begin
     if n_elements(guess) eq 0 and n_elements(init_img) gt 0 then begin
        w=where(hdunames eq init_img, count)
        if (count gt 1) then print, 'multiple init images found in input file, discarding all but first one.'
        if (count ge 1) then guess=mrdfits(input,hdunumbers[w[0]],header) ; TODO: use header coordinates etc.
     endif
  endif else begin
                                ; if init_img is a string, read it as
                                ; fits. Else check it is a 2d array
     sz = size(init_img) & nsz = N_ELEMENTS(sz) & typesz = sz[nsz-2]  
     if (typesz eq 7) then begin 
        guess=mrdfits(init_img,0,init_img_header)
        init_img=FILE_BASENAME(init_img) ; remove potentially harmful long dirname
     endif else begin 
        guess=init_img          ; a passed 2d array.
        delvar,init_img         ; avoid further us of init_img as if it was a filename.
     endelse
  end

; if guess is a cube, take 1st plane
if ((size(guess))[0] gt 2) then guess=guess[*,*,0,0,0,0,0]


; if rgl_prio is defined in header, then find this hdu, read image and use it as guess
  if ~dorgl_prio then begin
     prior=0                    ; a good starting point...
     if n_elements(rgl_prio) gt 0 then begin
        w=where(extname eq rgl_prio, count)
        if (count gt 1) then print, 'multiple init images found in input file, discarding.'
        if (count ge 1) then prior=mrdfits(input,w[0],header) ; TODO: use header coordinates etc.
     endif
  endif else prior=rgl_prio     ; defined in call 

;target: if not defined, take first one. If defined, find index, and
;error if target not found:
  if (target eq "*") then begin
     target=strtrim(oitarget[0].target,2)
     target_id=oitarget[0].target_id
     raep0=oitarget[0].raep0
     decep0=oitarget[0].decep0
  endif else begin
     targets=strtrim(oitarget.target,2)
     w=where(targets eq target, count)
     if (count lt 1) then message,"target not found in input file."
     target_id=oitarget[w[0]].target_id
     raep0=oitarget[w[0]].raep0
     decep0=oitarget[w[0]].decep0
  endelse
; warn (for the time being) that a complicated OIFITS is not pruned at
; output by the choice of only one object: 
if (n_elements(oitarget) gt 1) then message,/informational,"WARNING -- Output file will contain more HDUs than the selected target's ones."

; create table of correspondence: which wave correspond to vis2 and t3
  allinsnames=strarr(n_elements(oiwavearr))
  for i=0,n_elements(oiwavearr)-1 do allinsnames[i]=sxpar(*oiwaveheadarr[i],"INSNAME")

  t3inst=intarr(n_elements(oit3arr))
  vis2inst=intarr(n_elements(oivis2arr))

  for i=0,n_elements(oivis2arr)-1 do vis2inst[i]=(where(sxpar(*oivis2headarr[i],"INSNAME") eq allinsnames))[0]
  for i=0,n_elements(oit3arr)-1 do t3inst[i]=(where(sxpar(*oit3headarr[i],"INSNAME") eq allinsnames))[0]

; all values defined, read data
  data=wisard_oifits2data(input,targetname=target)

; there must be some data left to work with
  if (~finite(total(data.vis2err)) || ~finite(total(data.cloterr))) then message,"Only Flagged Data Available, Aborting."

  doWaveSubset=0
; simple wavelength selection. If selection is >= min max range, no
; subset is really asked for.
  if (wave_min lt wave_max) then begin ; which have probaly been read in input file
     w=where(data.wlen lt wave_min or data.wlen gt wave_max, count)
     if count gt 0 then begin ; there is indeed a subset asked!
        doWaveSubset=1
        w=where(data.wlen ge wave_min and data.wlen le wave_max, count)
        if count lt 1 then message,"no such wavelength range in data"
        data = data[w]
     endif
  end else begin                ; set values of wave_min and max to used values:
     wave_min = min(data.wlen)
     wave_max = max(data.wlen)
  end
  waveSubset=[wave_min,wave_max]

; update fov to be no greater than max fov given by telescope diameter
; and lambad_min
maxfov=(1.22*wave_min/(*oiarrayarr)[0].diameter)*180*3600.*1000./!DPI
if (fov gt maxfov or fov le 0) then begin
   print,'Setting FOV to (maximum) value of '+strtrim(string(maxfov, format='(F6.2)'),2)+'.'
   fov=maxfov
end 

; start reconstruction. TODO case wrt rgl_name/prior
  case regul of
     0: reconstruction = WISARD(data, NBITER = nbiter, threshold=threshold, GUESS = guess, FOV = fov*onemas, NP_MIN = np_min, AUX_OUTPUT = aux_output, /TOTVAR, oversampling=oversampling, positivity=positivity, display=display, USE_FLAGGED_DATA=use_flagged_data, simulated_data=issim) 
     
     1: begin
; if regul == 1 (psd) if no prior given, create the PSD of the documentation:
        if n_elements(prior) gt 1 then begin ; take centered ft or prior image as PSD
           psd=abs(shift(fft(prior),(size(prior))[1:2]/2+1))
        endif else begin
           distance = double(shift(dist(np_min),np_min/2,np_min/2))
           psd = 1D/((distance^3 > 1D) < 1D6)
        endelse
        reconstruction = WISARD(data, NBITER = nbiter, threshold=threshold, GUESS = guess, FOV = fov*onemas, NP_MIN = np_min, PSD=psd, AUX_OUTPUT = aux_output, oversampling=oversampling, positivity=positivity, display=display, USE_FLAGGED_DATA=use_flagged_data, simulated_data=issim)
     end

     2: begin
        scale = 1D/(NP_min)^2   ; factor for balance between regularization **** and in cost function
        delta = 1D
        reconstruction = WISARD(data, NBITER = nbiter, threshold=threshold, GUESS = guess, FOV = fov*onemas, NP_MIN = np_min, SCALE=scale, DELTA=delta, AUX_OUTPUT = aux_output, oversampling=oversampling, positivity=positivity, display=display, USE_FLAGGED_DATA=use_flagged_data, simulated_data=issim)
     end

     3: begin
        scale = 1D/(NP_min)^2   ; factor for balance between regularization **** and in cost function
        delta = 1D
        reconstruction = WISARD(data, NBITER = nbiter, threshold=threshold, GUESS = guess, FOV = fov*onemas, NP_MIN = np_min, SCALE=scale, DELTA=delta, /WHITE, AUX_OUTPUT = aux_output, oversampling=oversampling, positivity=positivity, display=display, USE_FLAGGED_DATA=use_flagged_data, simulated_data=issim)
     end

     4: begin
        if ~n_elements(fwhm) then fwhm=10                  ; pixels
        if ~n_elements(mu_support) then mu_support=10.0    ; why not?
        reconstruction = WISARD(data, NBITER = nbiter, threshold=threshold, GUESS = guess, FOV = fov*onemas, NP_MIN = np_min, MEAN_O=prior, MU_SUPPORT = mu_support, FWHM = fwhm, AUX_OUTPUT = aux_output, oversampling=oversampling, positivity=positivity, display=display, USE_FLAGGED_DATA=use_flagged_data, simulated_data=issim)
     end
     
     else: message,"regularization not yet supported, FIXME!"
  endcase

; prepare output HDU
  FXADDPAR,outhead,'XTENSION','BINTABLE'
  FXADDPAR,outhead,'EXTNAME','IMAGE-OI INPUT PARAM'
  FXADDPAR,outhead,'TARGET',target
  FXADDPAR,outhead,'WAVE_MIN',wave_min
  FXADDPAR,outhead,'WAVE_MAX',wave_max
  FXADDPAR,outhead,'USE_VIS', 'F'
  FXADDPAR,outhead,'USE_VIS2','T'
  FXADDPAR,outhead,'USE_T3', 'T'
  if (n_elements(init_img) gt 0) then FXADDPAR,outhead,'INIT_IMG',init_img ; the init image file passed
  FXADDPAR,outhead,'MAXITER',nbiter
  FXADDPAR,outhead,'RGL_NAME',regul_name[regul]
  if (n_elements(scale) gt 0 ) then FXADDPAR,outhead,'SCALE',scale
  if (n_elements(delta) gt 0 ) then FXADDPAR,outhead,'DELTA',delta
  FXADDPAR,outhead,'THRESHOL',threshold
  FXADDPAR,outhead,'NP_MIN',np_min
  FXADDPAR,outhead,'OVERSAMP',oversampling
  FXADDPAR,outhead,'FOV',fov,'Field of View (mas)'

; add correct values to images headers
  sz=size(aux_output.guess)
  nx=sz[1]
  ny=sz[2]
  if (n_elements(init_img) gt 0) then FXADDPAR,imagehead,'EXTNAME',init_img else FXADDPAR,imagehead,'EXTNAME','wisard_initimage'
  if (n_elements(init_img) gt 0) then FXADDPAR,imagehead,'HDUNAME',init_img else FXADDPAR,imagehead,'HDUNAME','wisard_initimage'
  FXADDPAR,imagehead,'CTYPE1','RA---SIN'
  FXADDPAR,imagehead,'CTYPE2','DEC--SIN'
  FXADDPAR,imagehead,'CRPIX1',nx/2
  FXADDPAR,imagehead,'CRPIX2',ny/2
  FXADDPAR,imagehead,'CROTA1',0.0
  FXADDPAR,imagehead,'CROTA2',0.0
  FXADDPAR,imagehead,'CRVAL1',raep0
  FXADDPAR,imagehead,'CRVAL2',decep0
  FXADDPAR,imagehead,'CDELT1',fov/3600./1000D/nx
  FXADDPAR,imagehead,'CDELT2',fov/3600./1000D/ny
  FXADDPAR,imagehead,'CUNIT1','deg'
  FXADDPAR,imagehead,'CUNIT2','deg'
  FXADDPAR,imagehead,'EQUINOX',2000.0
;main header: the reconstructed image
  sz=size(aux_output.x)
  nx=sz[1]
  ny=sz[2]
  FXADDPAR,main_header,'HDUNAME','wisard_image'
  FXADDPAR,main_header,'CTYPE1','RA---SIN'
  FXADDPAR,main_header,'CTYPE2','DEC--SIN'
  FXADDPAR,main_header,'CRPIX1',nx/2
  FXADDPAR,main_header,'CRPIX2',ny/2
  FXADDPAR,main_header,'CROTA1',0.0
  FXADDPAR,main_header,'CROTA2',0.0
  FXADDPAR,main_header,'CRVAL1',raep0
  FXADDPAR,main_header,'CRVAL2',decep0
  FXADDPAR,main_header,'CDELT1',fov/3600./1000D/nx
  FXADDPAR,main_header,'CDELT2',fov/3600./1000D/ny
  FXADDPAR,main_header,'CUNIT1','deg'
  FXADDPAR,main_header,'CUNIT2','deg'
  FXADDPAR,main_header,'EQUINOX',2000.0

; write output file. Put output image in main header!
  mwrfits,aux_output.x,output,main_header,/create,/silent,/no_copy,/no_comment
; create a dummy structure to make mwrfits write a binary extension
  dummystruct={DUMMY: 0.0d}
  mwrfits,dummystruct,output,outhead,/silent,/no_copy,/no_comment
  mwrfits,aux_output.guess,output,imagehead,/silent,/no_copy,/no_comment

  mwrfits,oitarget,output,targethead,/silent,/no_copy,/no_comment

; a good program would only write results related to current target! TODO!
  for i=0,n_elements(oiarrayarr)-1 do mwrfits,*(oiarrayarr[i]),output,*(oiarrayheadarr[i]),/silent,/no_copy,/no_comment

; do not write unknown arrays.
;  for i=0,n_elements(oiotherarr)-1 do mwrfits,*(oiotherarr[i]),output,*(oiotherheadarr[i]),/silent,/no_copy,/no_comment

  goodinsnamelist=''
; write structures updated with computed values: select only sections
; with target and pass corresponding wavelength for subset.
; memorize corresponding insnames to write only the corresponding
; wavelengths or wavelength subsets.
  for i=0,n_elements(oivis2arr)-1 do begin
     outhead=*oivis2headarr[i]
     if (doWaveSubset) then begin
        vis2subset=add_model_oivis2( *oivis2arr[i] , *oiwavearr[vis2inst[i]], aux_output, outhead, use_target=target_id, wsubs=waveSubset, operators=aux_output.operators)
     endif else begin
        vis2subset=add_model_oivis2( *oivis2arr[i] , *oiwavearr[vis2inst[i]], aux_output, outhead, use_target=target_id, operators=aux_output.operators)
     endelse
     if ( n_elements(vis2subset) gt 0 ) then begin
        ;add insname to 'good' insname list
        insname=strtrim(sxpar(outhead,"INSNAME"),2)
        goodinsnamelist=[goodinsnamelist,insname] ; add
        goodinsnamelist = goodinsnamelist[UNIQ(goodinsnamelist, SORT(goodinsnamelist))] ; sort
        col_of_flag=where(strtrim(tag_names(vis2subset),2) eq "FLAG", count) ; necessary for mwrfits.
        mwrfits,vis2subset, output, outhead,/silent,/no_copy,/no_comment,logical_cols=col_of_flag+1
     endif

  endfor

  for i=0,n_elements(oit3arr)-1 do begin
     outhead=*oit3headarr[i]
     if (doWaveSubset) then begin
        t3subset=add_model_oit3( *oit3arr[i], *oiwavearr[vis2inst[i]], aux_output, outhead, use_target=target_id, wsubs=waveSubset, operators=aux_output.operators )
     endif else begin 
        t3subset=add_model_oit3( *oit3arr[i], *oiwavearr[vis2inst[i]], aux_output, outhead, use_target=target_id , operators=aux_output.operators)
     endelse
    if ( n_elements(t3subset) gt 0 ) then begin ; no adding 'insname' as wisard uses both vis2 and t3. All relevant insnames have already been found.
        col_of_flag=where(strtrim(tag_names(t3subset),2) eq "FLAG", count)
        mwrfits,t3subset, output, outhead,/silent,/no_copy,/no_comment,logical_cols=col_of_flag+1
     endif
  endfor

  for i=0,n_elements(oiwavearr)-1 do begin ; write only relevant instrument tables
     insname=strtrim(sxpar(*(oiwaveheadarr[i]),"INSNAME"),2) ; avoid problems with blanks.
     w=where(goodinsnamelist eq insname, count)
     if (count gt 0) then begin
        if (doWaveSubset) then begin                                           ; write only wavelength subset for specific insname
           ww=where((*oiwavearr[i]).eff_wave ge waveSubset[0] and  (*oiwavearr[i]).eff_wave le waveSubset[1], count)
           if (count gt 0) then mwrfits,(*oiwavearr[i])[ww],output,*(oiwaveheadarr[i]),/silent,/no_comment
        endif else begin
           mwrfits,*(oiwavearr[i]),output,*(oiwaveheadarr[i]),/silent,/no_comment ; NO_COPY would invalidate the use of the table afterwards!
        endelse
     endif
  endfor



  if (~n_elements(interactive)) then begin
    print,'----------------ACKNOWLEDGEMENTS------------------------------'
    print,"If WISARD was helpful for your research, please add this sentence in the acknowledgement section of your articles:"
    print,"``This research has made use of the Jean-Marie Mariotti Center \textsc{WISARD} image reconstruction utility \footnote{Available at http://www.jmmc.fr/wisard}.''"
    print,"and cite the two following refereed papers in the body of your paper:"
    print,' (1) S. Meimon, L. M. Mugnier, and G. Le Besnerais, "Self-calibration approach for optical long-baseline interferometry imaging", J. Opt. Soc. Am. A, 26(1):108-120, 2009.'
    print,' (2) S. Meimon, L. M. Mugnier and G. Le Besnerais,``A convex approximation of the likelihood in optical interferometry'', J. Opt. Soc. Am. A (November 2005).'
     exit,status=0
  endif
end
