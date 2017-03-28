function add_model_oivis2,vis2,wave,aux_output,header, use_target=tid, wsubs=wsubs, operators=operators

; this functions adds the model-related columns to a OIVIS2 table,
; selecting only the target TOD and the wavelength subset wsubs. It
; returns a !NULL if the resulting table would be empty (no target or
; now corresponding wavelengths

  if (n_elements(tid) gt 0) then begin
     good=where(vis2.target_id eq tid, count)
     if count lt 1 then return,!NULL
     if count eq n_elements(vis2) then vis2=vis2[good]
  endif

  lambda=wave.eff_wave
  nwave=n_elements(lambda)
  ntimes=n_elements(vis2)

; safe way: discard eventual columns NS_MODEL_XXX, copy others and
; create our own NS_MODEL_XXX columns. Doing this, we must edit the
; header to remove TFORMx,TTYPEx,TUNITx values which will be updated
; by mwrfits.
  nhlines=n_elements(header)
  good=~strmatch(header,'TUNIT*')
  good=(good and ~strmatch(header,'TTYPE*'))
  good=(good and ~strmatch(header,'TFORM*'))
  w=where(good eq 1, count)
  if (count gt 0) then header=header[w] else header=''

  names=tag_names(vis2)
  ttag=where(names ne 'NS_MODEL_VIS2' and names ne 'NS_MODEL_VIS2ERR', count)
  ; silly way to create all other columns in the struct:
  temp_vis2=create_struct(names[ttag[0]],vis2[0].(ttag[0]))
  for itag=1,n_elements(ttag)-1 do temp_vis2=create_struct(temp_vis2,names[ttag[itag]],vis2[0].(ttag[itag])) ; everything minus above columns.
  ; header is input header minus TFORMx,TTYPEx,TUNITx 


  if (nwave gt 1) then new_vis2=create_struct('NS_MODEL_VIS2',fltarr(nwave)) else new_vis2=create_struct('NS_MODEL_VIS2',0.0)
  new_vis2=replicate(create_struct(temp_vis2,new_vis2),n_elements(vis2))

  for itag=0,n_elements(ttag)-1 do new_vis2.(itag)=vis2.(ttag[itag]) ; populate temp_vis2's first columns of new_vis2.

;;;; Computation of spatial frequencies
;freqs_u = operators._B#aux_output.freqs_u
;freqs_v = operators._B#aux_output.freqs_v
; compute model:
; replicate ucoord to nwave, replicate lambda to ntimes, flatten,
; divide for spatial freqs:
  uvector=transpose(cmreplicate(vis2.ucoord,nwave))
  vvector=transpose(cmreplicate(vis2.vcoord,nwave))

  lambdavector=cmreplicate(lambda,ntimes)
  freqs_u=reform(uvector/lambdavector,nwave*ntimes)
  freqs_v=reform(vvector/lambdavector,nwave*ntimes)

  np=(size(aux_output.x))[1]

  H=WISARD_MAKE_H(FREQS_U=freqs_u, FREQS_V=freqs_v,$
                  FOV = aux_output.fov, NP_MIN = -np,$
                  NP_OUTPUT = NPOUT, STEP_OUTPUT = step_output)
; NP may not be np_min, although returned image (x) is only np_min x
; np_min. congrid is needed. (but probably a bad idea!)
;  congrid does not preserve flux: renormalize!
  norm_x = reform(aux_output.x,np*np)
  norm_x /= total(norm_x)

  achix = reform(H#norm_x) ;
  new_vis2.ns_model_vis2=reform(reform(abs2(achix),nwave,ntimes))

  if (n_elements(wsubs) eq 0) then begin
     ; add relevant units to header
     names=tag_names(new_vis2)
     for itag=1,n_elements(names)-1 do begin
        name=names[itag]
        tunitstr='TUNIT'+strtrim(itag+1,2)
        if (name eq 'TIME') then   FXADDPAR,header,tunitstr,'s'
        if (name eq 'MJD') then   FXADDPAR,header,tunitstr,'day'
        if (name eq 'INT_TIME') then   FXADDPAR,header,tunitstr,'s'
        if (name eq 'UCOORD') then   FXADDPAR,header,tunitstr,'m'
        if (name eq 'VCOORD') then   FXADDPAR,header,tunitstr,'m'
     end
     return, new_vis2
  endif
; or: apply wavelength subset if necessary

  wl=where(lambda ge wsubs[0] and lambda le wsubs[1], count)
  if (count le 0) then return,!NULL
  
  lambda=lambda(wl)
  nwave=count

  ; must create a new table/structure since spectral dimension has changed:
  names=tag_names(new_vis2)
  ttag=where(names ne "VIS2DATA" and names ne "VIS2ERR" and names ne "NS_MODEL_VIS2" and names ne "FLAG", count)
  ; silly way to create all other columns in the struct:
  vis2subset=create_struct(names[ttag[0]],new_vis2[0].(ttag[0]))
  for itag=1,n_elements(ttag)-1 do vis2subset=create_struct(vis2subset,names[ttag[itag]],new_vis2[0].(ttag[itag]))

  ; now add wavelength subsets:
  if (nwave gt 1) then vis2addsubset=create_struct('VIS2DATA', dblarr(nwave), 'VIS2ERR', dblarr(nwave), 'NS_MODEL_VIS2',fltarr(nwave), 'FLAG', bytarr(nwave)) else vis2addsubset=create_struct('VIS2DATA', 0.0d, 'VIS2ERR', 0.0d, 'NS_MODEL_VIS2',0.0,'FLAG',OB)
   new_vis2subset=replicate(create_struct(vis2subset,vis2addsubset),n_elements(new_vis2))
   ; populate new_vis2subset for non-wavelength tags
   for itag=0,n_elements(ttag)-1 do new_vis2subset.(itag)=new_vis2.(ttag[itag])
   ; populate for wavelength subsets
   new_vis2subset.vis2data=new_vis2.vis2data[wl]
   new_vis2subset.vis2err=new_vis2.vis2err[wl]
   new_vis2subset.ns_model_vis2=new_vis2.ns_model_vis2[wl]
   new_vis2subset.flag=new_vis2.flag[wl]
   ; add relevant units to header
   names=tag_names(new_vis2subset)
   for itag=1,n_elements(names)-1 do begin
      name=names[itag]
      tunitstr='TUNIT'+strtrim(itag+1,2)
      if (name eq 'TIME') then   FXADDPAR,header,tunitstr,'s'
      if (name eq 'MJD') then   FXADDPAR,header,tunitstr,'day'
      if (name eq 'INT_TIME') then   FXADDPAR,header,tunitstr,'s'
      if (name eq 'UCOORD') then   FXADDPAR,header,tunitstr,'m'
      if (name eq 'VCOORD') then   FXADDPAR,header,tunitstr,'m'
   end
   return, new_vis2subset
end
