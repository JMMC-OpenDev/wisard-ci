;
;+
;
; NAME: MODEL2OIFITS
;
; PURPOSE: 
;  creates an OImaging OIFITS (OIFITS with special entries added)
;  using a regular OIFITS file as input and a fits file containing a
;  reconstructed image plus a special table MODEL-VISIBILITIES written
;  by a image reconstruction prgram that do not want to write
;  complicated oifits by itself. Hopefully, but not tested by this
;  procedure, the reconstructed image and the u,v spatial frequencies
;  listed in the MODEL-VISIBILITIES tabel are those of the input OIFITS
;  file for the selected target (which defaults to the first target
;  found)
;       
; CATEGORY:
; Optical Image Reconstruction Helper, OIFITS.
; CALLING SEQUENCE: 
;  MODEL2OIFITS,INPUT_OIFITS,MODEL_FITS[, OUTPUT_OIFITS],TARGET=TARGET
; KEYWORD PARAMETERS:
; TARGET: pass a target name (in case there are more than one target
; in the oifits file).
;
; LICENCE:
; Copyright (C) 2018, G. Duvert unless explicitly mentioned in files.
;
; This program is free software; you can redistribute it and/or modify  
; it under the terms of the GNU General Public License as published by  
; the Free Software Foundation; either version 2 of the License, or     
; (at your option) any later version.                                   
;
;-

pro fill_in_model_vis2,mvis2,vis2,wave,model
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
         if (count gt 0) then mvis2[ivis].ns_model_vis2[iwave]=(smodel[ww[0]].model_visamp)^2
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

pro fill_in_model_t3,mt3,t3,wave,model
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

; ----------------------- THE MAIN ROUTINE ------------------------------------------

pro model2oifits,input,model,output,target=target

  ; if /display option, we are interactive
  @ "wisard_common.pro"
  wisard_is_interactive =  keyword_set(display)
  @ "wisard_catch_noniteractive.pro"

  dotarget=n_elements(target) ne 0

;; defaults in absence of passed values. Normally these should be read in one of the headers.
  if ~dotarget    then target="*"

;; FITS_INFO is robust in case of strange tables (no column table). 
  fits_info,input, n_ext=next, extname=extname, /silent
  if (next le 4) then message,"Input file is probably not an OIFITS file at all!"
  extname=strtrim(extname,2)
;;create list of hdu numbers; and list of hdunames(if any) for further use:
  hdunames=replicate('',next)
  hdunumbers=indgen(next)

;read first header
  hdu0=mrdfits(input,0,main_header)
;examine others
  for i=1,next do begin
     if extname[i] eq "OI_T3" then begin
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

; read simple fits from imagin program: image (primary),
; and  model-visibilities

  fits_info,model, n_ext=mext, extname=m_extname, /silent
  if (mext lt 1) then message,"Model file is probably not what we want!"
  m_extname=strtrim(m_extname,2)
;;create list of hdu numbers; and list of hdunames(if any) for further use:
  m_hdunames=replicate('',mext)
  m_hdunumbers=indgen(mext)

;read first header
  image=mrdfits(model,0,model_main_header)
;examine others, get visibilities and input/output params.
  for i=1,mext do begin
     ;print,"Model HDU: ",i, " extname: ",m_extname[i]
     if (m_extname[i] eq "IMAGE-OI MODEL VISIBILITIES" or m_extname[i] eq "MODEL-VISIBILITIES") then begin ; found in some old version of mira
        mvis = mrdfits(model,i,m_header)
     endif else if (m_extname[i] eq  "IMAGE-OI INPUT PARAM") then begin ; only one input param.
        inputparam = mrdfits(model,i,inputparamhead)
     endif else if (m_extname[i] eq  "IMAGE-OI OUTPUT PARAM") then begin ; only one output param.
        outputparam = mrdfits(model,i,outputparamhead)
     endif
  endfor
; if no "model", trouble:
  if (n_elements(mvis) lt 1) then message,"Wrong model file, exiting"


; write output file. Put output image in main header!
;main header: the reconstructed image. Just a copy of main header of 'model':
  sz=size(image)
  nx=sz[1]
  ny=sz[2]
  mwrfits,image,output,model_main_header,/create,/silent,/no_copy,/no_comment

; write params (header). trick is that we NEED to pass a structure.
  if (n_elements(inputparam)) then mwrfits,{dummyin:0.0},output,inputparamhead,/silent,/no_copy,/no_comment
  if (n_elements(outputparam)) then mwrfits,{dummyout:0.0},output,outputparamhead,/silent,/no_copy,/no_comment

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
      vis2subset=insert_model_oivis2( *oivis2arr[i] , *oiwavearr[vis2inst[i]], mvis, outhead, use_target=target_id)
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
      t3subset=insert_model_oit3( *oit3arr[i], *oiwavearr[t3inst[i]], mvis, outhead, use_target=target_id)
      if ( n_elements(t3subset) gt 0 ) then begin ; no adding 'insname' as wisard uses both vis2 and t3. All relevant insnames have already been found.
         col_of_flag=where(strtrim(tag_names(t3subset),2) eq "FLAG", count)
         mwrfits,t3subset, output, outhead,/silent,/no_copy,/no_comment,logical_cols=col_of_flag+1
      endif
   endfor

  for i=0,n_elements(oiwavearr)-1 do begin ; write only relevant instrument tables
     insname=strtrim(sxpar(*(oiwaveheadarr[i]),"INSNAME"),2) ; avoid problems with blanks.
     w=where(goodinsnamelist eq insname, count)
     if (count gt 0) then begin
       mwrfits,*(oiwavearr[i]),output,*(oiwaveheadarr[i]),/silent,/no_comment ; NO_COPY would invalidate the use of the table afterwards!
     endif
  endfor

exit,status=0
end

