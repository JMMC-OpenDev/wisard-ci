function add_model_oit3,t3,wave,aux_output,header, use_target=tid, wsubs=wsubs, operators=operators

; this functions adds the model-related columns to a OIVIS2 table,
; selecting only the target TOD and the wavelength subset wsubs. It
; returns a !NULL if the resulting table would be empty (no target or
; now corresponding wavelengths
 
  if (n_elements(tid) gt 0) then begin
     good=where(t3.target_id eq tid, count)
     if count lt 1 then return,!NULL
     if count ne n_elements(t3) then t3=t3[good] ; select only good!
  endif 


  lambda=wave.eff_wave
  nwave=n_elements(lambda)
  ntimes=n_elements(t3)

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

  names=tag_names(t3)
  ttag=where(names ne 'NS_MODEL_T3AMP' and names ne 'NS_MODEL_T3AMPERR' and names ne 'NS_MODEL_T3PHI' and names ne 'NS_MODEL_T3PHIERR' , count)
  ; silly way to create all other columns in the struct:
  temp_t3=create_struct(names[ttag[0]],t3[0].(ttag[0])) 
  for itag=1,n_elements(ttag)-1 do temp_t3=create_struct(temp_t3,names[ttag[itag]],t3[0].(ttag[itag])) ; everything minus above columns.

  if (nwave gt 1) then new_t3=create_struct('NS_MODEL_T3AMP',fltarr(nwave),'NS_MODEL_T3PHI',fltarr(nwave)) else new_t3=create_struct('NS_MODEL_T3AMP',0.0,'NS_MODEL_T3PHI',0.0)
  new_t3=replicate(create_struct(temp_t3,new_t3),n_elements(t3))

  for itag=0,n_elements(ttag)-1 do new_t3.(itag)=t3.(ttag[itag]) ; populate 'temp_t3' first columns of new_t3.

; compute model:
; replicate ucoord to nwave, replicate lambda to ntimes, flatten,
; divide for spatial freqs:
  u1vector=transpose(cmreplicate(t3.u1coord,nwave)) ; u coordinate of baseline AB of the triangle (m)
  u2vector=transpose(cmreplicate(t3.u2coord,nwave)) ; u coordinate of baseline BC of the triangle (m)
  v1vector=transpose(cmreplicate(t3.v1coord,nwave)) ; v coordinate of baseline AB of the triangle (m)
  v2vector=transpose(cmreplicate(t3.v2coord,nwave)) ; v coordinate of baseline BC of the triangle (m)

  lambdavector=cmreplicate(lambda,ntimes)

  freqs_u1=reform(u1vector/lambdavector,nwave*ntimes)
  freqs_v1=reform(v1vector/lambdavector,nwave*ntimes)
  freqs_u2=reform(u2vector/lambdavector,nwave*ntimes)
  freqs_v2=reform(v2vector/lambdavector,nwave*ntimes)

; for all 3 visibilities
; convention : if NP_min<0 it is a desired NP (forced) 
; NP may not be np_min, although returned image (x) is only np_min x
; np_min. congrid is needed. (but probably a bad idea!)
;  congrid does not preserve flux: renormalize!

  np=(size(aux_output.x))[1]

  H=WISARD_MAKE_H(FREQS_U=freqs_u1, FREQS_V=freqs_v1,$
                  FOV = aux_output.fov, NP_MIN = -np,$
                  NP_OUTPUT = NPOUT, STEP_OUTPUT = step_output)
  norm_x = reform(aux_output.x,np*np)
  norm_x /= total(norm_x)
  achix1 = reform(H#norm_x) ;

  H=WISARD_MAKE_H(FREQS_U=freqs_u2, FREQS_V=freqs_v2,$
                  FOV = aux_output.fov, NP_MIN = -np,$
                  NP_OUTPUT = NPOUT, STEP_OUTPUT = step_output)

  norm_x = reform(aux_output.x,np*np)
  norm_x /= total(norm_x)
  achix2 = reform(H#norm_x) ;

  H=WISARD_MAKE_H(FREQS_U=(freqs_u1+freqs_u2), FREQS_V=(freqs_v1+freqs_v2),$ ; baseline AC
                  FOV = aux_output.fov, NP_MIN = -np,$
                  NP_OUTPUT = NPOUT, STEP_OUTPUT = step_output)

  norm_x = reform(aux_output.x,np*np)
  norm_x /= total(norm_x)
  achix3 = reform(H#norm_x)

  tripleproduct=achix1*achix2*conj(achix3)

  new_t3.ns_model_t3amp=reform(reform(abs(tripleproduct),nwave,ntimes)) 
  new_t3.ns_model_t3phi=reform(reform(atan(tripleproduct,/phase),nwave,ntimes))
  ;convert phases to degrees!
  new_t3.ns_model_t3phi*=180.0d/!DPI

  if (n_elements(wsubs) eq 0) then begin
     ; add relevant units to header
     names=tag_names(new_t3)
     for itag=1,n_elements(names)-1 do begin
        name=names[itag]
        tunitstr='TUNIT'+strtrim(itag+1,2)
        if (name eq 'TIME') then   FXADDPAR,header,tunitstr,'s'
        if (name eq 'MJD') then   FXADDPAR,header,tunitstr,'day'
        if (name eq 'INT_TIME') then   FXADDPAR,header,tunitstr,'s'
        if (name eq 'U1COORD') then   FXADDPAR,header,tunitstr,'m'
        if (name eq 'V1COORD') then   FXADDPAR,header,tunitstr,'m'
        if (name eq 'U2COORD') then   FXADDPAR,header,tunitstr,'m'
        if (name eq 'V2COORD') then   FXADDPAR,header,tunitstr,'m'
        if (name eq 'T3PHI') then   FXADDPAR,header,tunitstr,'deg'
        if (name eq 'T3PHIERR') then   FXADDPAR,header,tunitstr,'deg'
        if (name eq 'NS_MODEL_T3PHI') then   FXADDPAR,header,tunitstr,'deg'
     end
     return, new_t3
  endif
; or: apply wavelength subset if necessary

  wl=where(lambda ge wsubs[0] and lambda le wsubs[1], count)
  if (count le 0) then return,!NULL

  lambda=lambda(wl)
  nwave=count

  ; must create a new table/structure since spectral dimension has changed:
  names=tag_names(new_t3)

  ttag=where(names ne 'NS_MODEL_T3AMP' and names ne 'NS_MODEL_T3PHI' and names ne 'T3AMP' and names ne 'T3AMPERR' and names ne 'T3PHI' and names ne 'T3PHIERR' and names ne 'FLAG', count)
  ; silly way to create all other columns in the struct:
  t3subset=create_struct(names[ttag[0]],new_t3[0].(ttag[0]))
  for itag=1,n_elements(ttag)-1 do t3subset=create_struct(t3subset,names[ttag[itag]],new_t3[0].(ttag[itag]))

  ; now add wavelength subsets:
  if (nwave gt 1) then t3addsubset=create_struct('T3AMP', dblarr(nwave), 'T3AMPERR', dblarr(nwave), 'T3PHI',  dblarr(nwave), 'T3PHIERR', dblarr(nwave),'NS_MODEL_T3AMP', fltarr(nwave), 'NS_MODEL_T3PHI',  fltarr(nwave),  'FLAG', bytarr(nwave)) else t3addsubset=create_struct('T3AMP',  0.0d, 'T3AMPERR', 0.0d, 'T3PHI',  0.0d, 'T3PHIERR', 0.0d,'NS_MODEL_T3AMP', 0.0, 'NS_MODEL_T3PHI',  0.0, 'FLAG',OB)
   new_t3subset=replicate(create_struct(t3subset,t3addsubset),n_elements(new_t3))
   ; populate new_t3subset for non-wavelength tags
   for itag=0,n_elements(ttag)-1 do new_t3subset.(itag)=new_t3.(ttag[itag])
   ; populate for wavelength subsets
   new_t3subset.t3amp=new_t3.t3amp[wl]
   new_t3subset.t3amperr=new_t3.t3amperr[wl]
   new_t3subset.t3phi=new_t3.t3phi[wl]
   new_t3subset.t3phierr=new_t3.t3phierr[wl]
   new_t3subset.ns_model_t3amp=new_t3.ns_model_t3amp[wl]
   new_t3subset.ns_model_t3phi=new_t3.ns_model_t3phi[wl]
   new_t3subset.flag=new_t3.flag[wl]
   ; add relevant units to header
   names=tag_names(new_t3subset)
   for itag=1,n_elements(names)-1 do begin
      name=names[itag]
      tunitstr='TUNIT'+strtrim(itag+1,2)
      if (name eq 'TIME') then   FXADDPAR,header,tunitstr,'s'
      if (name eq 'MJD') then   FXADDPAR,header,tunitstr,'day'
      if (name eq 'INT_TIME') then   FXADDPAR,header,tunitstr,'s'
      if (name eq 'U1COORD') then   FXADDPAR,header,tunitstr,'m'
      if (name eq 'V1COORD') then   FXADDPAR,header,tunitstr,'m'
      if (name eq 'U2COORD') then   FXADDPAR,header,tunitstr,'m'
      if (name eq 'V2COORD') then   FXADDPAR,header,tunitstr,'m'
      if (name eq 'T3PHI') then   FXADDPAR,header,tunitstr,'deg'
      if (name eq 'T3PHIERR') then   FXADDPAR,header,tunitstr,'deg'
      if (name eq 'NS_MODEL_T3PHI') then   FXADDPAR,header,tunitstr,'deg'
   end
   return, new_t3subset
end

