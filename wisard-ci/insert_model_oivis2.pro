pro  fill_in_model_vis2,mvis2,vis2,wave,model
                                ; checks based on dimension of mvis2:
  sz=size(mvis2.ns_model_vis2)
  nvis=(sz[0] gt 1)?sz[2]:sz[1]
  nwave=(sz[0] gt 1)?sz[1]:1
  if (sz[-1] ne nwave*nvis) then message,"internal error in fill_in_model_vis2"
  for ivis=0,nvis-1 do begin    ; ivis same order between mvis2 and vis2.vis2data
     ucoord=vis2[ivis].ucoord
     vcoord=vis2[ivis].vcoord
     ; mira uses only ucoord >O so we convert in -u -v if u < 0
     if (ucoord lt 0.0) then begin ucoord=-ucoord & vcoord=-vcoord & end
     w=where( abs(model.ucoord-ucoord) lt 1d-3 and abs(model.vcoord-vcoord lt 1d-3), count)
     if (count gt 0) then begin
       smodel=model[w] 
       for iwave=0,nwave-1 do begin ; iwave same order with wave.eff_wave
         wlen=wave[iwave].eff_wave*1D9 ; mira in nm.
         ww=where( abs(smodel.eff_wave-wlen) lt 1e-2, count)
         if (count gt 0) then mvis2[ivis].ns_model_vis2[iwave]=smodel[ww[0]].model_visamp
       endfor
     endif
  endfor
end

function insert_model_oivis2,vis2,wave,mvis,header, use_target=tid

; this functions adds the model-related columns to a OIVIS2 table,
; selecting only the target TID and the wavelength subset wsubs. It
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

  fill_in_model_vis2,new_vis2,vis2,wave,mvis

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
end
