function add_model_oivis2,vis2,wave,data, use_target=tid, wsubs=wsubs

  if (n_elements(tid) gt 0) then begin
     good=where(vis2.target_id eq tid, count)
     if count lt 1 then return,vis2
     if count eq n_elements(vis2) then begin
        goodvis2=vis2
        all = 1
     endif else begin
        goodvis2=vis2[good]
        all = 0
     endelse
  endif else begin
     goodvis2=vis2
     all = 1
  end


  addvis2=1B
  lambda=wave.eff_wave
  nwave=n_elements(lambda)
  ntimes=n_elements(vis2)
  ngoodtimes=n_elements(goodvis2)

  names=tag_names(vis2)
  w=where(names eq "NS_MODEL_VIS2", count)
  if (count gt 0) then addvis2=0B 
  if (addvis2) then begin
     if (nwave gt 1) then new_vis2=create_struct('NS_MODEL_VIS2',fltarr(nwave),'NS_MODEL_VIS2ERR',fltarr(nwave)) else new_vis2=create_struct('NS_MODEL_VIS2',0.0,'NS_MODEL_VIS2ERR',0.0)
     new_vis2=replicate(create_struct(vis2[0],new_vis2),n_elements(vis2))
     ; populate new_vis2 
     for itag=0,n_tags(vis2)-1 do new_vis2.(itag)=vis2.(itag)
  endif else new_vis2=vis2

; replicate ucoord to nwave, replicate lambda to ngoodtimes, flatten,
; divide for spatial freqs:
  uvector=transpose(cmreplicate(goodvis2.ucoord,nwave))
  vvector=transpose(cmreplicate(goodvis2.vcoord,nwave))
  lambdavector=cmreplicate(lambda,ngoodtimes)
  freqs_u=reform(uvector/lambdavector,nwave*ngoodtimes)
  freqs_v=reform(vvector/lambdavector,nwave*ngoodtimes)
  H=WISARD_MAKE_H(FREQS_U=freqs_u, FREQS_V=freqs_v,$
                  FOV = data.fov, NP_MIN = data.np_min,$
                  NP_OUTPUT = NP, STEP_OUTPUT = step_output)
; NP may not be np_min, although returned image (x) is only np_min x
; np_min. congrid is needed. (but probably a bad idea!)
;  congrid does not preserve flux: renormalize!
  norm_x = reform(congrid(data.x,np,np),np*np)
  norm_x /= total(norm_x)

  achix = reform(H#norm_x) ;
  if (all) then new_vis2.ns_model_vis2=reform(reform(abs2(achix),nwave,ngoodtimes)) else new_vis2.ns_model_vis2[good]=reform(reform(abs2(achix),nwave,ngoodtimes))

  if (n_elements(wsubs) eq 0) then return, new_vis2
; or: apply wavelength subset if necessary

  wl=where(lambda ge wsubs[0] and lambda le wsubs[1], count)
  if (count le 0) then begin 
     message,/informational,"ERROR in wavelength subset values in function add_model_oivis2()."
     exit,status=1
  endif
  lambda=lambda(wl)
  nwave=count
  ; must create a new table/structure
  names=tag_names(new_vis2)
  ttag=where(names ne "VIS2DATA" and names ne "VIS2ERR" and names ne "NS_MODEL_VIS2" and names ne "NS_MODEL_VIS2ERR" and names ne "FLAG", count)
  ; silly way to create all other columns in the struct:
  vis2subset=create_struct(names[ttag[0]],new_vis2[0].(ttag[0]))
  for itag=1,n_elements(ttag)-1 do vis2subset=create_struct(vis2subset,names[ttag[itag]],new_vis2[0].(ttag[itag]))

  ; now add wavelength subsets:
  if (nwave gt 1) then vis2addsubset=create_struct('VIS2DATA', dblarr(nwave), 'VIS2ERR', dblarr(nwave), 'NS_MODEL_VIS2',fltarr(nwave),'NS_MODEL_VIS2ERR',fltarr(nwave), 'FLAG', bytarr(nwave)) else vis2addsubset=create_struct('VIS2DATA', 0.0d, 'VIS2ERR', 0.0d, 'NS_MODEL_VIS2',0.0,'NS_MODEL_VIS2ERR',0.0,'FLAG',OB)
   new_vis2subset=replicate(create_struct(vis2subset,vis2addsubset),n_elements(new_vis2))
   ; populate new_vis2subset for non-wavelength tags
   for itag=0,n_elements(ttag)-1 do new_vis2subset.(itag)=new_vis2.(ttag[itag])
   ; populate for wavelength subsets
   new_vis2subset.vis2data=new_vis2.vis2data[wl]
   new_vis2subset.vis2err=new_vis2.vis2err[wl]
   new_vis2subset.ns_model_vis2=new_vis2.ns_model_vis2[wl]
   new_vis2subset.ns_model_vis2err=new_vis2.ns_model_vis2err[wl]
   new_vis2subset.flag=new_vis2.flag[wl]
   return, new_vis2subset
end
