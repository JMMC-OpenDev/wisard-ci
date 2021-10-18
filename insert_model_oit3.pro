pro  fill_in_model_t3,mt3,t3,wave,model
                                ; checks based on dimension of mvis2:
  sz=size(mt3.ns_model_t3amp)
  nt3=(sz[0] gt 1)?sz[2]:sz[1]
  nwave=(sz[0] gt 1)?sz[1]:1
  if (sz[-1] ne nwave*nt3) then message,"internal error in fill_in_model_t3"
  for it3=0,nt3-1 do begin    ; it3 same order between mt3 and t3.t3amp/t3phi
     u1coord=t3[it3].u1coord
     v1coord=t3[it3].v1coord
     ; mira uses only ucoord >O so we convert in -u -v if u < 0
     if (u1coord lt 0.0) then begin u1coord=-u1coord & v1coord=-v1coord & end
     u2coord=t3[it3].u2coord
     v2coord=t3[it3].v2coord
     ; mira uses only ucoord >O so we convert in -u -v if u < 0
     if (u2coord lt 0.0) then begin u2coord=-u2coord & v2coord=-v2coord & end
     u3coord=t3[it3].u2coord+t3[it3].u1coord
     v3coord=t3[it3].v2coord+t3[it3].v1coord
     if (u3coord lt 0.0) then begin u3coord=-u3coord & v3coord=-v3coord & end

     w1=where( abs(model.ucoord-u1coord) lt 1d-3 and abs(model.vcoord-v1coord lt 1d-3), count1)
     w2=where( abs(model.ucoord-u2coord) lt 1d-3 and abs(model.vcoord-v2coord lt 1d-3), count2)
     w3=where( abs(model.ucoord-u3coord) lt 1d-3 and abs(model.vcoord-v3coord lt 1d-3), count3)
     if (count1 gt 0 and count2 gt 0 and count3 gt 0 ) then begin
       smodel1=model[w1] 
       smodel2=model[w2] 
       smodel3=model[w3] 
       for iwave=0,nwave-1 do begin ; iwave same order with wave.eff_wave
         wlen=wave[iwave].eff_wave*1D9 ; mira in nm.
         ww1=where( abs(smodel1.eff_wave-wlen) lt 1e-2, count1)
         ww2=where( abs(smodel2.eff_wave-wlen) lt 1e-2, count2)
         ww3=where( abs(smodel3.eff_wave-wlen) lt 1e-2, count3)
         if (count1 gt 0 and count2 gt 0 and count3 gt 0) then begin
           mt3[it3].ns_model_t3amp[iwave]=smodel1[ww1[0]].model_visamp*smodel2[ww2[0]].model_visamp*smodel3[ww3[0]].model_visamp
           mt3[it3].ns_model_t3phi[iwave]=smodel1[ww1[0]].model_visphi+smodel2[ww2[0]].model_visphi-smodel3[ww3[0]].model_visphi
         endif
       endfor
     endif
  endfor
end

function insert_model_oit3,t3,wave,mvis,header, use_target=tid

; this functions adds the model-related columns to a OIT3 table,
; selecting only the target TID and the wavelength subset wsubs. It
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

  fill_in_model_t3,new_t3,t3,wave,mvis

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
end
