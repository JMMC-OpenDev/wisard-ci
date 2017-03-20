function add_model_oit3,t3,wave,data,use_target=tid, wsubs=wsubs

  if (n_elements(tid) gt 0) then begin
     good=where(t3.target_id eq tid, count)
     if count lt 1 then return,t3
     if count eq n_elements(t3) then begin
        goodt3=t3
        all = 1
     endif else begin
        goodt3=t3[good]
        all = 0
     endelse
  endif else begin
     goodt3=t3
     all = 1
  end


  addt3=1B
  lambda=wave.eff_wave
  nwave=n_elements(lambda)
  ntimes=n_elements(t3)
  ngoodtimes=n_elements(goodt3)
  

  names=tag_names(t3)
  w=where(names eq "NS_MODEL_T3AMP", count)
  if (count gt 0) then addt3=0B 
  if (addt3) then begin
     if (nwave gt 1) then new_t3=create_struct('NS_MODEL_T3AMP',fltarr(nwave),'NS_MODEL_T3AMPERR',fltarr(nwave),'NS_MODEL_T3PHI',fltarr(nwave),'NS_MODEL_T3PHIERR',fltarr(nwave)) else new_t3=create_struct('NS_MODEL_T3AMP',0.0,'NS_MODEL_T3AMPERR',0.0,'NS_MODEL_T3PHI',0.0,'NS_MODEL_T3PHIERR',0.0)
     new_t3=replicate(create_struct(t3[0],new_t3),n_elements(t3))
     ; populate new_t3 
     for itag=0,n_tags(t3)-1 do new_t3.(itag)=t3.(itag)
  endif else new_t3=t3


; replicate ucoord to nwave, replicate lambda to ngoodtimes, flatten,
; divide for spatial freqs:
  u1vector=transpose(cmreplicate(goodt3.u1coord,nwave))
  u2vector=transpose(cmreplicate(goodt3.u2coord,nwave))
  v1vector=transpose(cmreplicate(goodt3.v1coord,nwave))
  v2vector=transpose(cmreplicate(goodt3.v2coord,nwave))
  lambdavector=cmreplicate(lambda,ngoodtimes)
  freqs_u1=reform(u1vector/lambdavector,nwave*ngoodtimes)
  freqs_v1=reform(v1vector/lambdavector,nwave*ngoodtimes)
  freqs_u2=reform(u2vector/lambdavector,nwave*ngoodtimes)
  freqs_v2=reform(v2vector/lambdavector,nwave*ngoodtimes)

; for all 3 visibilities
; convention : si NP_min<0 alors c'est le NP voulu exactement 
; NP may not be np_min, although returned image (x) is only np_min x
; np_min. congrid is needed. (but probably a bad idea!)
;  congrid does not preserve flux: renormalize!

  H=WISARD_MAKE_H(FREQS_U=freqs_u1, FREQS_V=freqs_v1,$
                  FOV = data.fov, NP_MIN = data.np_min,$
                  NP_OUTPUT = NP, STEP_OUTPUT = step_output)
  norm_x = reform(congrid(data.x,np,np),np*np)
  norm_x /= total(norm_x)
  achix1 = reform(H#norm_x) ;

  H=WISARD_MAKE_H(FREQS_U=freqs_u2, FREQS_V=freqs_v2,$
                  FOV = data.fov, NP_MIN = data.np_min,$
                  NP_OUTPUT = NP, STEP_OUTPUT = step_output)
  norm_x = reform(congrid(data.x,np,np),np*np)
  norm_x /= total(norm_x)
  achix2 = reform(H#norm_x) ;

  H=WISARD_MAKE_H(FREQS_U=(freqs_u2+freqs_u1), FREQS_V=(freqs_v2+freqs_v1),$
                  FOV = data.fov, NP_MIN = data.np_min,$
                  NP_OUTPUT = NP, STEP_OUTPUT = step_output)
  norm_x = reform(congrid(data.x,np,np),np*np)
  norm_x /= total(norm_x)
  achix3 = reform(H#norm_x)

  tripleproduct=achix1*achix2*conj(achix3)
  if (all) then new_t3.ns_model_t3amp=reform(reform(abs(tripleproduct),nwave,ngoodtimes)) else new_t3.ns_model_t3amp[good]=reform(reform(abs(tripleproduct),nwave,ngoodtimes))

  if (all) then new_t3.ns_model_t3phi=reform(reform(atan(tripleproduct,/phase),nwave,ngoodtimes)) else new_t3.ns_model_t3phi[good]=reform(reform(atan(tripleproduct,/phase),nwave,ngoodtimes))
  ;convert phases to degrees!
  if (all) then new_t3.ns_model_t3phi*=180.0d/!DPI else new_t3.ns_model_t3phi[good]*=180.0d/!DPI

  if (n_elements(wsubs) eq 0) then return, new_t3
; or: apply wavelength subset if necessary

  wl=where(lambda ge wsubs[0] and lambda le wsubs[1], count)
  if (count le 0) then begin 
     message,/informational,"ERROR in wavelength subset values in function add_model_oivis2()."
     exit,status=1
  endif
  lambda=lambda(wl)
  nwave=count

  ; must create a new table/structure
  names=tag_names(new_t3)

  ttag=where(names ne 'NS_MODEL_T3AMP' and names ne 'NS_MODEL_T3AMPERR' and names ne 'NS_MODEL_T3PHI' and names ne 'NS_MODEL_T3PHIERR' and names ne 'T3AMP' and names ne 'T3AMPERR' and names ne 'T3PHI' and names ne 'T3PHIERR' and names ne 'FLAG', count)
  ; silly way to create all other columns in the struct:
  t3subset=create_struct(names[ttag[0]],new_t3[0].(ttag[0]))
  for itag=1,n_elements(ttag)-1 do t3subset=create_struct(t3subset,names[ttag[itag]],new_t3[0].(ttag[itag]))

  ; now add wavelength subsets:
  if (nwave gt 1) then t3addsubset=create_struct('NS_MODEL_T3AMP', fltarr(nwave), 'NS_MODEL_T3AMPERR', fltarr(nwave), 'NS_MODEL_T3PHI',  fltarr(nwave), 'NS_MODEL_T3PHIERR', fltarr(nwave), 'T3AMP',  dblarr(nwave), 'T3AMPERR', dblarr(nwave), 'T3PHI',  dblarr(nwave), 'T3PHIERR', dblarr(nwave), 'FLAG', bytarr(nwave)) else t3addsubset=create_struct('NS_MODEL_T3AMP', 0.0, 'NS_MODEL_T3AMPERR', 0.0, 'NS_MODEL_T3PHI',  0.0, 'NS_MODEL_T3PHIERR', 0.0, 'T3AMP',  0.0d, 'T3AMPERR', 0.0d, 'T3PHI',  0.0d, 'T3PHIERR', 0.0d,'FLAG',OB)
   new_t3subset=replicate(create_struct(t3subset,t3addsubset),n_elements(new_t3))
   ; populate new_t3subset for non-wavelength tags
   for itag=0,n_elements(ttag)-1 do new_t3subset.(itag)=new_t3.(ttag[itag])
   ; populate for wavelength subsets
   new_t3subset.t3amp=new_t3.t3amp[wl]
   new_t3subset.t3amperr=new_t3.t3amperr[wl]
   new_t3subset.t3phi=new_t3.t3phi[wl]
   new_t3subset.t3phierr=new_t3.t3phierr[wl]
   new_t3subset.ns_model_t3amp=new_t3.ns_model_t3amp[wl]
   new_t3subset.ns_model_t3amperr=new_t3.ns_model_t3amperr[wl]
   new_t3subset.ns_model_t3phi=new_t3.ns_model_t3phi[wl]
   new_t3subset.ns_model_t3phierr=new_t3.ns_model_t3phierr[wl]
   new_t3subset.flag=new_t3.flag[wl]
   return, new_t3subset
end

