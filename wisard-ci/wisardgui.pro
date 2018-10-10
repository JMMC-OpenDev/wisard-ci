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
;     SIMULATED_DATA=SIMULATED_DATA,
;     USE_FLAGGED_DATA=USE_FLAGGED_DATA,
;     RECONSTRUCTED_IMAGE=MYIMAGE,SCALE=SCALE,DELTA=DELTA,/HELP
;
; KEYWORD PARAMETERS:
;    INPUT: input OIFITS or (better) OImaging OIFITS file (that contains
;    already guess image and parameters). 
;    Otherwise, use Keywords to enter these parameters
;
;    OUTPUT: product: an OImaging OIFITS file name. Image is in main header.
;
;    TARGET: the object name (in case there are many in the input file)
;    NBITER: max number of iterations (50 by default)
;    NP_MIN: minimum number of resels to reconstruct. default computed
;    internally 
;    REGUL = one of ['TOTVAR','L1L2','L1L2WHITE','PSD','SOFT_SUPPORT']
;    INIT_IMG: Guess start image (fits). Supposedly not mandatory
;    according to the doc although unescapable in practice.
;
;    WAVERANGE=[min,max] to limit reconstruction to this wave
;    interval. Values in micrometers. Must include at least ONE
;    channel.
;    RECONSTRUCTED_IMAGE=xx the image will also be returned in the
;    variable xx. Interactive mode only. Can be used as an imput
;    image, as init_img reads also IDL/GDL variables. 
;
;    The (other) Keywords are best described in the documentation, see:
;    WISARD:
;    http://www.mariotti.fr/doc/approved/JMMC-MAN-2500-0001.pdf
;    OImaging interface:
;    http://www.mariotti.fr/doc/approved/JRA4Fp7Report2016.pdf
;
; EXAMPLE:
;    wisardgui,"inputdata/2004/2004-FKV1137.fits","output_of_example.fits", $
;    nbiter=1000,fov=16,np_min=64,waverange=[1.0,5.0],reconstructed_image=myimage, $
;    regul='L1L2WHITE',/display,/simulated_data
;
; NOTE: 
;    At the moment, the model visibility will be forced in ALL wavelengths, not only the
;    ones that were used (case where WAVERANGE was specified). FIXME!
; LICENCE:
; Copyright (C) 2018, G. Duvert unless explicitly mentioned in files.
;
; This program is free software; you can redistribute it and/or modify  
; it under the terms of the GNU General Public License as published by  
; the Free Software Foundation; either version 2 of the License, or     
; (at your option) any later version.                                   
;
;-

pro wisardgui,input,output,target=target,threshold=threshold,nbiter=nbiter,fov=fov,np_min=np_min,regul=regul,positivity=positivity,oversampling=oversampling,init_img=passed_init_img,rgl_prio=rgl_prio,display=display,mu_support=mu_support, fwhm=fwhm, waverange=waverange, simulated_data=issim, use_flagged_data=use_flagged_data,reconstructed_img=reconstructed_img,scale=scale,delta=delta,_extra=ex,help=help

; if /display option, we are interactive
@ "wisard_common.pro"
wisard_is_interactive =  keyword_set(display)
@ "wisard_catch_noniteractive.pro"

if keyword_set(help) then begin
 doc_library,"wisardgui"
 if (~(wisard_is_interactive)) then exit else return
end

  forward_function WISARD_OIFITS2DATA,WISARD ; GDL does not need that, IDL insists on it!

  wisard_ci_version=getenv('WISARD_CI_VERSION')
  if (strlen( wisard_ci_version) lt 1) then wisard_ci_version='Unversioned'
  if wisard_is_interactive then extra_comment="" else extra_comment=", batch mode"
  message,/informational,"Welcome to Wisard, version "+wisard_ci_version+", you have accepted the copyrights."
  message,/informational,"Wisard is running with "+strcompress(!CPU.TPOOL_NTHREADS)+" threads"+extra_comment+"."

  if wisard_is_interactive then device, decomposed=0, retain=2

  ; memorize passed line values
  dotarget=n_elements(target) ne 0 
  dothreshold=n_elements(threshold) ne 0
  donbiter=n_elements(nbiter) ne 0
  dofov=n_elements(fov) ne 0
  donp_min=n_elements(np_min) ne 0
  doregul=n_elements(regul) ne 0
  doinit_img=n_elements(passed_init_img) ne 0
  dorgl_prio=n_elements(rgl_prio) ne 0
  dowaverange=n_elements(waverange) ne 0
  doScale=n_elements(scale) ne 0
  doDelta=n_elements(delta) ne 0
  if dotarget then passed_target=target 
  if dothreshold then passed_threshold=threshold
  if donbiter then passed_nbiter=nbiter
  if dofov then passed_fov=fov
  if donp_min then passed_np_min = np_min
  if doregul then passed_regul=regul
  if dorgl_prio then passed_rgl_prio=rgl_prio
  if dowaverange then passed_waverange=waverange
  if doScale then passed_scale=scale
  if doDelta then passed_delta=delta



 ; after catch to handle this first error message if necessary. 
  if (~strmatch(!PATH,'*wisardlib*')) then begin
     b=routine_info('WISARDGUI',/source)
     c=file_dirname(b.path)
     path_to_Wisard=c+"/lib"
     !path=expand_path('+'+path_to_Wisard)+":"+!path
     if (~strmatch(!PATH,'*wisardlib*')) then message,"Unable to set PATH to wisard libraries. Please set up yourself"
  endif

;; define constants
  regul_name=['L1L2','L1L2WHITE','PSD','SOFT_SUPPORT','TOTVAR'] & default_regul=0 ; L1L2
  defsysv, '!DI', dcomplex(0, 1), 1 ;Define i; i^2=-1 as readonly system-variable
  onemas=1d-3*(!DPi/180D)/3600D     ; 1 mas in radian

;; non-OImaging options:
  doreconstructed_img=arg_present(reconstructed_img)
  doovers=n_elements(oversampling) ne 0
  doposit=n_elements(positivity) ne 0
  if ~doovers     then oversampling=1 else oversampling = fix (oversampling) ; oversampling
  if ~doposit     then positivity=1   else positivity = fix (positivity) ; quoted: "misplaced curiosity"
  
;; defaults in absence of passed values. Will be updated by the ones
;; in the input file, which are superseded by the eventual arguments passed.
  target='*'
  threshold=1d-6
  nbiter=50
  fov=14.0D
  np_min=32L
  regul=default_regul
  wave_min = -1
  wave_max = -1
  guess = 0
  init_img = ''
  prior = 0
  find_init_img = 0
  find_prior_img = 0

;; First read input file and get start image & parameters that overread defaults.
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
        fits_read,input,dummy,input_params_header,/header_only,exten_no=i
; parse relevant parameters

; TARGET
        req_target=strtrim(sxpar(input_params_header,"TARGET"),2) ; avoid problems with blanks.
        if req_target ne '0' then target=req_target

; WAVE_MIN
        req_wave_min=sxpar(input_params_header,"WAVE_MIN")
        if req_wave_min ne '0' then wave_min=req_wave_min

; WAVE_MAX
        req_wave_max=sxpar(input_params_header,"WAVE_MAX")
        if req_wave_max ne '0' then wave_max=req_wave_max

; USE_VIS,USE_VIS2,USE_T3: ignored

; INIT_IMG
        req_init_img=strtrim(sxpar(input_params_header,"INIT_IMG"),2) ; avoid problems with blanks.
        if req_init_img ne '0' then begin
           init_img=req_init_img ; it is the name of a HDUNAME.
           find_init_img = 1 ; will force to find it
           message,/inform,"... OIMaging file requests "+init_img+" as init image"
        endif

; MAXITER
        req_nbiter=sxpar(input_params_header,"MAXITER")
        if req_nbiter gt 0 then nbiter=req_nbiter

; RGL_NAME
        req_rgl_name=strtrim(sxpar(input_params_header,"RGL_NAME"),2) ; avoid problems with blanks.
        if req_rgl_name ne '0' then begin
           w=where(strtrim(req_rgl_name,2) eq regul_name, count)
           if count lt 1 then regul=default_regul else regul=w[0]
        endif

; RGL_WGHT,RGL_ALPHA,RGL_BETA: ignored at the moment

; RGL_PRIO
        req_rgl_prio=strtrim(sxpar(input_params_header,"RGL_PRIO"),2) ; avoid problems with blanks.
        if req_rgl_prio ne '0' then begin
           rgl_prio=req_rgl_prio ; it is the name of a HDUNAME.
           find_prior_img=1
           message,/inform,"... OIMaging file requests "+rgl_prio+" as regularization prior"
        endif
; AUTO_WGHT ?
; NP_MIN
        req_np_min=sxpar(input_params_header,"NP_MIN") 
        if req_np_min gt 0 then np_min=req_np_min

; FLUXERR -> threshold
        req_fluxerr=sxpar(input_params_header,"FLUXERR") 
        if req_fluxerr gt 0 then threshold=double(req_fluxerr)/(np_min)^2

; FOV
        req_fov=sxpar(input_params_header,"FOV") 
        if req_fov gt 0 then fov=req_fov
        
; WISARD: SCALE: WARNING: This is not the same logic. As SCALE default
; is best computed from the image and NP, wee need to know if a SCALE
; was forced either by the commandline or the OIMaging file.
        if (~doScale) then begin   ; factor for balance between regularization **** and in cost function
           req_scale=sxpar(input_params_header,"SCALE") 
           if req_scale gt 0 then begin ; meaning <= 0 will provide the default.
              scale=req_scale
              doScale = 1
           endif
        endif
; WISARD: DELTA: WARNIG: Same logic as for SCALE/
        if (~doDelta) then begin   ; factor for balance between regularization **** and in cost function
           req_delta=sxpar(input_params_header,"DELTA") 
           if req_delta gt 0 then begin; meaning <= 0 will provide the default.
              delta=req_delta
              doDelta = 1
           endif
        endif

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

; if init_img is defined in header, then find this hdu, read image and
; use it as guess. remember HDU as we  will copy it.
  if find_init_img then begin
     w=where(hdunames eq init_img, count)
     if (count gt 1) then print, 'multiple init images found in input file, discarding all but first one.'
     if (count ge 1) then begin 
        hdu_init_img=hdunumbers[w[0]]
        guess=mrdfits(input,hdu_init_img,init_img_header) ; TODO: use header coordinates etc.
        message,/inform,'... found init image "'+init_img+'" at HDU #'+strtrim(string(hdu_init_img),2)
     endif else find_init_img = 0
  endif

; if rgl_prio is defined in header, then find this hdu, read image and use it as guess
  if find_prior_img then begin
     w=where(extname eq rgl_prio, count)
     if (count gt 1) then print, 'multiple prior images found in input file, discarding.'
     if (count ge 1) then begin
        hdu_prior_img=hdunumbers[w[0]]
        prior=mrdfits(input,hdu_prior_img,prior_img_header) ; TODO: use header coordinates etc.
        message,/inform,'... found prior image "'+rgl_prio+'" at HDU #'+strtrim(string(hdu_prior_img),2)
     endif else find_prior_img = 0
  endif

;; input file read complete. Overwrite values with passed values if any:

  if doinit_img then begin
     find_init_img = 0 ; since it's going to be replaced
     ; if passed_init_img is a string, read it as fits. Else check it is a 2d array
     sz = size(passed_init_img) & nsz = N_ELEMENTS(sz) & typesz = sz[nsz-2]  
     if (typesz eq 7) then begin 
        guess=mrdfits(passed_init_img,0,init_img_header)
        init_img=FILE_BASENAME(passed_init_img) ; remove potentially harmful long dirname
        message,/inform,'... using external fits file "'+init_img+'" as init image' 
     endif else begin 
        guess=passed_init_img   ; a passed 2d array.
        init_img="internal_image"
        message,/inform,"... using passed array as init image" 
     endelse
  endif

  if dorgl_prio then begin
     find_prior_impg = 0
     ; if prior is a string, read it as fits. Else check it is a 2d array
     sz = size(passed_rgl_prio) & nsz = N_ELEMENTS(sz) & typesz = sz[nsz-2]  
     if (typesz eq 7) then begin 
        prior=mrdfits(passed_rgl_prio,0,rgl_prio_header)
        rgl_prio=FILE_BASENAME(rgl_prio) ; remove potentially harmful long dirname
        message,/inform,'... using external fits file "'+rgl_prio+'" as prior' 
     endif else begin 
        prior=passed_rgl_prio          ; a passed 2d array.
        rgl_prio="internal_image"
        message,/inform,"... using passed array as prior" 
     endelse
  endif


  if dotarget    then target=passed_target
  if dothreshold then threshold=double(passed_threshold); convergence threshold
  if donbiter    then nbiter=fix(passed_nbiter) ; number of iterations.
  if dofov       then fov = double(passed_fov) ; Field of View of reconstructed image (14*14 marcsec^2 here)
  if donp_min    then np_min = fix(passed_np_min) ; width of reconstructed image in pixels (same along x and y axis).
  if doregul     then begin ; L1L2 by default
      w=where(strtrim(strupcase(regul_name),2) eq passed_regul, count) ; regul is a name here
      if count lt 1 then regul=default_regul else regul=w[0]    ; regul is now an integer: 0:totvar, 1:psd, 2:l1l2, 3:l1l2_white, 4:support 
  endif

  if (dowaverange) then begin
     type = SIZE(waverange, /TYPE)
     if (type eq 7) then result=execute('waverange='+waverange) ; to convert string to array
     wave_min=waverange[0]*1D-6
     wave_max=waverange[1]*1D-6
  endif

  message,/informational,"Using Regularization: "+regul_name[regul]

; if guess is a cube, take 1st plane
  if ((size(guess))[0] gt 2) then guess=guess[*,*,0,0,0,0,0]


;target: if still not defined, take first one. If defined, find index, and
;error if target not found:
  if (target eq "*") then begin
     target=strtrim(oitarget[0].target,2)
     target_id=oitarget[0].target_id
     raep0=oitarget[0].raep0
     decep0=oitarget[0].decep0
  endif else begin
     targets=strtrim(oitarget.target,2)
     w=where(targets eq target, count)
     if (count lt 1) then begin
        message,/inform,'ERROR: target "'+target+'" not found in input file.'
        if (~(wisard_is_interactive)) then exit else return
     endif
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

; all values defined,  get new data structure
  masterDataArray=wisard_oifits2data(input,targetname=target,_extra=ex)
; drop void data in masterData 
  www=where(ptr_valid(masterDataArray),numberOfConfigurations)
  masterDataArray=masterDataArray[www]
; individually check an trim each data structure
  someDataIsAvailable=0
  someRangeIsAvailable=(wave_min eq wave_max) ; 1 if we do not select any wave range.
; total data wavelength range: compute, intializing to some good
; value:
  total_wave_min=min((*masterDataArray[0]).wlen)
  total_wave_max=max((*masterDataArray[0]).wlen)
; use last (greatest number of baselines) not to break compatibility
  for idata=0,numberOfConfigurations-1 do begin
     data=*masterDataArray[idata]
     
; there must be some data left to work with
     if ( total(data.vis2err,/nan) ne 0 && total(data.cloterr,/nan) ne 0 ) then someDataIsAvailable++
     
; simple wavelength selection. If selection is >= min max range, no
; subset is really asked for. wave_min==wave_mask==-1 if no range has
; been asked for.
     if (wave_min lt wave_max) then begin ; which have probaly been read in input file
        w=where(data.wlen lt wave_min or data.wlen gt wave_max, count)
        if count gt 0 then begin ; there is indeed a subset asked!
           w=where(data.wlen ge wave_min and data.wlen le wave_max, count)
           if count gt 0 then begin
              data = data[w]
              masterDataArray[idata]=PTR_NEW(data)
              someRangeIsAvailable++
           endif
        endif
     endif else begin
        total_wave_min=min([total_wave_min,data.wlen])        
        total_wave_max=max([total_wave_max,data.wlen])        
     endelse
  endfor


  if ( ~someDataIsAvailable) then message,"ERROR: Only Flagged Data Available, Aborting."
  if ( ~someRangeIsAvailable) then message,"ERROR: no data found within requested wavelength range, Aborting. (request is probably smaller than a single instrument channel?)"

; if wave_min,max is still defaulting to -1, use total wavelength values:
  if (wave_min eq wave_max) then begin
     wave_min=total_wave_min
     wave_max=total_wave_max
  end

;TBC AS THERE MAY BE MORE THAN 1 ARRAY NOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
; update fov to be no greater than max fov given by telescope diameter
; and lambad_min
maxfov=(1.22*wave_min/((*(oiarrayarr[0])).diameter)[0])*180*3600.*1000./!DPI
if (fov gt maxfov or fov le 0) then begin
   print,'Setting FOV to (maximum) value of '+strtrim(string(maxfov, format='(F6.2)'),2)+'.'
   fov=maxfov
end 

; format some help/debug line
 commandline='Line equivalent of your command: wisardgui,/display,"'+input+'","'+output+'",target="'+target+'",threshold='+strtrim(string(threshold),2)+',nbiter='+strtrim(string(nbiter),2)+',fov='+strtrim(string(fov),2)+',np_min='+strtrim(string(np_min),2)+',regul="'+regul_name[regul]+'",waverange=['+strtrim(string(wave_min*1E6),2)+','+strtrim(string(wave_max*1E6),2)+']'
 if (doinit_img) then commandline+='",init_img="'+init_img
 if (dorgl_prio) then commandline+='",rgl_prio="'+rgl_prio+'"'
 if (doScale) then commandline+=',scale='+strtrim(string(scale),2)
 if (doDelta) then commandline+=',delta='+strtrim(string(delta),2)
  print,commandline

; start reconstruction. TODO case wrt rgl_name/prior
; TOTVAR is just delta=very_small. 
  case regul of

     0: begin ; L1L2
        if (~doScale) then scale = 1D/(NP_min)^2   ; factor for balance between regularization **** and in cost function
        if (~doDelta) then delta = 1D
        reconstruction = WISARD(masterDataArray, NBITER = nbiter, threshold=threshold, GUESS = guess, FOV = fov*onemas, NP_MIN = np_min, SCALE=scale, DELTA=delta, AUX_OUTPUT = aux_output, oversampling=oversampling, positivity=positivity, display=display, USE_FLAGGED_DATA=use_flagged_data, simulated_data=issim, _extra=ex)
     end

     1: begin ; L1L2WHITE
        if (~doScale) then scale = 1D/(NP_min)^2   ; factor for balance between regularization **** and in cost function
        if (~doDelta) then delta = 1D
        reconstruction = WISARD(masterDataArray, NBITER = nbiter, threshold=threshold, GUESS = guess, FOV = fov*onemas, NP_MIN = np_min, SCALE=scale, DELTA=delta, /WHITE, AUX_OUTPUT = aux_output, oversampling=oversampling, positivity=positivity, display=display, USE_FLAGGED_DATA=use_flagged_data, simulated_data=issim, _extra=ex)
     end

     2: begin ; PSD
; if regul == 1 (psd) if no prior given, create the PSD of the documentation:
        if n_elements(prior) gt 1 then begin ; take centered ft or prior image as PSD
           psd=abs(shift(fft(prior),(size(prior))[1:2]/2+1))
        endif else begin
           distance = double(shift(dist(np_min),np_min/2,np_min/2))
           psd = 1D/((distance^3 > 1D) < 1D6)
        endelse
        reconstruction = WISARD(masterDataArray, NBITER = nbiter, threshold=threshold, GUESS = guess, FOV = fov*onemas, NP_MIN = np_min, PSD=psd, AUX_OUTPUT = aux_output, oversampling=oversampling, positivity=positivity, display=display, USE_FLAGGED_DATA=use_flagged_data, simulated_data=issim, _extra=ex)
     end

     3: begin ; SOFT_SUPPORT
        if ~n_elements(fwhm) then fwhm=10                  ; pixels
        if ~n_elements(mu_support) then mu_support=10.0    ; why not?
        reconstruction = WISARD(masterDataArray, NBITER = nbiter, threshold=threshold, GUESS = guess, FOV = fov*onemas, NP_MIN = np_min, MEAN_O=prior, MU_SUPPORT = mu_support, FWHM = fwhm, AUX_OUTPUT = aux_output, oversampling=oversampling, positivity=positivity, display=display, USE_FLAGGED_DATA=use_flagged_data, simulated_data=issim, _extra=ex)
     end

     4: begin ; TOTVAR
        if (~doScale) then scale = 1D/(NP_min)^2   ; factor for balance between regularization **** and in cost function
        delta = scale/1D9 ; TOTVAR is L1L2 with delta infinitesimal.
        reconstruction = WISARD(masterDataArray, NBITER = nbiter, threshold=threshold, GUESS = guess, FOV = fov*onemas, NP_MIN = np_min, SCALE=scale, DELTA=delta, AUX_OUTPUT = aux_output, oversampling=oversampling, positivity=positivity, display=display, USE_FLAGGED_DATA=use_flagged_data, simulated_data=issim, _extra=ex)
     end
     
     else: message,"ERROR: regularization not yet supported, FIXME!"
  endcase

; main header: the reconstructed image. always present.
  sz=size(aux_output.x)
  nx=sz[1]
  ny=sz[2]
  snp=strtrim(string(np_min),2)
  sit=strtrim(string(aux_output.iter),2)

  reconstructed_image_hduname=target+'_'+snp+'x'+snp+'_'+sit
  FXADDPAR,main_header,'HDUNAME',reconstructed_image_hduname
  FXADDPAR,main_header,'CTYPE1','RA---TAN'
  FXADDPAR,main_header,'CTYPE2','DEC--TAN'
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
; initialize output file. Put output image in main header!
  mwrfits,aux_output.x,output,main_header,/create,/silent,/no_copy,/no_comment

; output parameters. always present.
  FXADDPAR,ouput_params_header,'XTENSION','BINTABLE'
  FXADDPAR,ouput_params_header,'EXTNAME','IMAGE-OI OUTPUT PARAM'
  FXADDPAR,ouput_params_header,'TARGET',target
  FXADDPAR,ouput_params_header,'WAVE_MIN',wave_min
  FXADDPAR,ouput_params_header,'WAVE_MAX',wave_max
  FXADDPAR,ouput_params_header,'LAST_IMG',reconstructed_image_hduname
  FXADDPAR,ouput_params_header,'MAXITER',nbiter
  FXADDPAR,ouput_params_header,'RGL_NAME',regul_name[regul]
  if doScale then FXADDPAR,ouput_params_header,'SCALE',scale
  if doDelta then FXADDPAR,ouput_params_header,'DELTA',delta
  FXADDPAR,ouput_params_header,'THRESHOL',threshold
  FXADDPAR,ouput_params_header,'NP_MIN',np_min
  FXADDPAR,ouput_params_header,'FOV',fov,'Field of View (mas)'
  FXADDPAR,ouput_params_header,'SOFTWARE','WISARD','IR software name'
  FXADDPAR,ouput_params_header,'VERSION',wisard_ci_version,'version of software'

  dummystruct={DUMMY: 0.0d}; we need a dummy structure to make mwrfits write a binary extension
  mwrfits,dummystruct,output,ouput_params_header,/silent,/no_copy,/no_comment

; input params only if there was an image inside. In which case we
; copy the initial image too, under its initial name
  if find_init_img then begin
     mwrfits,{aaa:0.0},output,input_params_header,/silent,/no_copy,/no_comment

; reinterpret input parameters and add correct values to images headers
     
     sz=size(aux_output.guess)
     nx=sz[1]
     ny=sz[2]
     FXADDPAR,init_image_header,'EXTNAME','INITIAL_IMAGE'
     FXADDPAR,init_image_header,'HDUNAME',init_img
     FXADDPAR,init_image_header,'CTYPE1','RA---TAN'
     FXADDPAR,init_image_header,'CTYPE2','DEC--TAN'
     FXADDPAR,init_image_header,'CRPIX1',nx/2
     FXADDPAR,init_image_header,'CRPIX2',ny/2
     FXADDPAR,init_image_header,'CROTA1',0.0
     FXADDPAR,init_image_header,'CROTA2',0.0
     FXADDPAR,init_image_header,'CRVAL1',raep0
     FXADDPAR,init_image_header,'CRVAL2',decep0
     FXADDPAR,init_image_header,'CDELT1',fov/3600./1000D/nx
     FXADDPAR,init_image_header,'CDELT2',fov/3600./1000D/ny
     FXADDPAR,init_image_header,'CUNIT1','deg'
     FXADDPAR,init_image_header,'CUNIT2','deg'
     FXADDPAR,init_image_header,'EQUINOX',2000.0
     mwrfits,aux_output.guess,output,init_image_header,/silent,/no_copy,/no_comment

  endif

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


     doWaveSubset=0
; note: the doWaveSubset trick is not correctly implemented in the
; subroutine and is temporarily ignored if fact. Meaning that the
; model visibility will be forced in ALL wavelengths, not only the
; ones that were used.

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
; note: the doWaveSubset trick is not correctly implemented in the subroutine
     if (doWaveSubset) then begin
        t3subset=add_model_oit3( *oit3arr[i], *oiwavearr[t3inst[i]], aux_output, outhead, use_target=target_id, wsubs=waveSubset, operators=aux_output.operators )
     endif else begin 
        t3subset=add_model_oit3( *oit3arr[i], *oiwavearr[t3inst[i]], aux_output, outhead, use_target=target_id , operators=aux_output.operators)
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
        ; note: the doWaveSubset trick is not correctly implemented in the subroutine
        if (doWaveSubset) then begin                                           ; write only wavelength subset for specific insname
           ww=where((*oiwavearr[i]).eff_wave ge waveSubset[0] and  (*oiwavearr[i]).eff_wave le waveSubset[1], count)
           if (count gt 0) then mwrfits,(*oiwavearr[i])[ww],output,*(oiwaveheadarr[i]),/silent,/no_comment
        endif else begin
           mwrfits,*(oiwavearr[i]),output,*(oiwaveheadarr[i]),/silent,/no_comment ; NO_COPY would invalidate the use of the table afterwards!
        endelse
     endif
  endfor

  if (~(wisard_is_interactive)) then begin
    print,'----------------ACKNOWLEDGEMENTS------------------------------'
    print,"If WISARD was helpful for your research, please add this sentence in the acknowledgement section of your articles:"
    print,"``This research has made use of the Jean-Marie Mariotti Center \textsc{WISARD} image reconstruction utility \footnote{Available at http://www.jmmc.fr/wisard}.''"
    print,"and cite the two following refereed papers in the body of your paper:"
    print,' (1) S. Meimon, L. M. Mugnier, and G. Le Besnerais, "Self-calibration approach for optical long-baseline interferometry imaging", J. Opt. Soc. Am. A, 26(1):108-120, 2009.'
    print,' (2) S. Meimon, L. M. Mugnier and G. Le Besnerais,``A convex approximation of the likelihood in optical interferometry'', J. Opt. Soc. Am. A (November 2005).'
     exit,status=0
  endif else if doreconstructed_img then reconstructed_img=aux_output.x 
end

