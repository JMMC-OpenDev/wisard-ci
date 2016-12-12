;+
; NAME:
;     wisard_build_data
; PURPOSE:
;     to sort non-friendly structure output by test.pro in a complete,
;     ordered, wisard-friendly structure
;
; CALLING SEQUENCE:
;     out=wisard_build_data(in)
;
; INPUTS:
;     aa
;
; OPTIONAL INPUTS:
;     VERBOSE=[0..2]   be more.. eh.. verbose as number increases 
;     /NOFLATTEN: keep spectral dimension as 3rd dimension in output
;     structure.
;     /USEVISDATA If data contains VISDATA complex values
;     (e.g. ASPRO simulations or AMBER data), keep them in output
;     struct, so as to have wisard_data2mdata use them instead of its
;     normal behaviour.
;
; OUTPUTS:
;
; RESTRICTIONS:
;       None.
;
; PROCEDURE:
; given a data struct containing all complete closures/visibilities at
; different times, for each time determine the number of telescopes
; involved ntel, order the closures in 123 124 125 234 235 345, ie,
; all the triplets with tel1, all the remaining with tel2, etc. keep
; all v2 associated.
;
; EXAMPLE:
;
; PROCEDURES USED:
;
; MODIFICATION HISTORY:
;     v0.0 2013    G. Duvert
;-

function wisard_build_data,data,verbose=verbose, noflatten=noflatten, usevisdata=usevisdata, numberoftels=numberoftels, first=first
;ON_ERROR,2
  if (n_elements(verbose) eq 0) then verbose=0
  nofirst=0 ; 'first' was given
  if (n_elements(first) eq 0) then begin nofirst=1 & first=0 & endif
  free_numberoftels=0
  if (n_elements(numberoftels) eq 0) then begin numberoftels=-1 & free_numberoftels=1 & endif
   warnHasIncompatibleConfigurations=0
   warnHasUncompleteConfigurationSet=0
   ntimes=(size(data,/dim))[0]  ; ntimes contain eventually nwavelength in p backends also if flatten (default)
   nelementsOfStructure=(size(*data[0],/dim))[1]
   hasOiVisData=(nelementsOfStructure eq 17 and keyword_set(usevisdata)) ;
   hasOiVis=(nelementsOfStructure ge 15)                                 ; to be updated at need

   ;initial dataset value: count nclots:
   initialsize=0 & for iTime=0L,ntimes-1 do initialsize+=(size(*data[iTime],/dim))[0]
   desiredSize=-1               ; define, but to a silly value
   strindex=0
   for iTime=0L,ntimes-1 do begin
      bundle=*data[iTime]
      nclot=(size(bundle,/dim))[0]
      tellist=intarr(3*nclot)
      for iClot=0L,nClot-1 do tellist[iClot*3:iClot*3+2]=*bundle[iClot,1]
      aa=histogram(tellist,binsize=1) ;
      ntels=n_elements(where(aa gt 0))
      ;; test and validate ntels passed by numberoftels argument
      if (numberoftels eq -1) then begin
         numberoftels=ntels
         message,/informational,"Step 2: Data looks like it is from a "+strcompress(string(numberoftels))+" telescopes array. If wrong change value with numberoftels=xx option."
      endif else begin
         if (ntels lt numberoftels) then begin 
            if (verbose gt 0) then print,format='(%"Data for time %f has %d telescopes, request of %d tels is impossible to fulfill")',*bundle[0,0],ntels,numberoftels
            goto, NOT_GOOD_ENOUGH ;
         endif else begin
            if (verbose gt 0 and ntels gt numberoftels) then print,format='(%"Data for time %f has %d telescopes, more than the %d requested")',*bundle[0,0],ntels,numberoftels
            ntels=numberoftels
         endelse
      endelse 
      nbases=ntels*(ntels-1)/2
      nclostheoric=ntels*(ntels-1)*(ntels-2)/6
      nindeptclostheoric=(ntels-1)*(ntels-2)/2
      ncoords=ntels-1

      if (verbose gt 0) then begin
         if (nClot lt nindeptclostheoric) then print,format='(%"Data for time %f has less closures (%d) than needed (%d), dropped.")',*bundle[0,0],nClot,nindeptclostheoric
         if (nClot gt nindeptclostheoric) then print,format='(%"Data for time %f has more closures (%d) than needed (%d).")',*bundle[0,0],nClot,nindeptclostheoric
         if (nClot eq nindeptclostheoric) then print,format='(%"Data for time %f has no redundant closures. Noise optimisation will not be possible.")',*bundle[0,0]
      endif

      ;; if we precise the reference telescope it means we want to check all configurations, so we need all possible closures. drop if not available
      if (nClot lt nclostheoric and nofirst eq 0 ) then goto, NOT_GOOD_ENOUGH

      mytels=tellist[uniq(tellist,SORT(tellist))] 
      ;;create an ordered triplet list. 
      mytels=shift(temporary(mytels),-first) 
      myclos=intarr(3,nindeptclostheoric)
      index=0
      mytel1=mytels[0]
      for J=1L,ntels-1 do begin 
         mytel2=mytels[j] 
         for K=J+1,ntels-1 do begin 
            mytel3=mytels[K] 
            myclos[*,index]=[mytel1,mytel2,mytel3] 
            index++ 
         endfor 
      endfor
      ;;start with actual closures at time iTime
      doneclot=intarr(nindeptclostheoric)
      for J=0L,nindeptclostheoric-1 do begin
         for iClot=0L,nclot-1 do begin
            the_triplet=*bundle[iClot,1]
            sign= oifits_triplet_compare(myclos[*,J],the_triplet)
            if (sign ne 0) then begin
               (*bundle[iClot,1])[0]=myclos[0,J]
               (*bundle[iClot,1])[1]=myclos[1,J]
               (*bundle[iClot,1])[2]=myclos[2,J]
               
               thet3phi=sign*[*(*bundle[iClot,2])]*!DPI/180D
               thet3phierr=[*(*bundle[iClot,3])]*!DPI/180D
               thet3flag=*(*bundle[iClot,4])
               www=where(thet3flag eq 70, count) & if (count GT 0) then thet3flag[www]=0
               www=where(thet3flag eq 84, count) & if (count GT 0) then thet3flag[www]=1
               
               wherenan=where(~finite(thet3phi), nancount) & if(nancount GT 0) then thet3flag[wherenan]=1  
               wherenan=where(~finite(thet3phierr), nancount) & if(nancount GT 0) then thet3flag[wherenan]=1  
               clot=(total(size(clot)) LT 1)?[thet3phi]:[[clot],[thet3phi]]
               cloterr=(total(size(cloterr)) LT 1)?[thet3phierr]:[[cloterr],[thet3phierr]]
               clotflag=(total(size(clotflag)) LT 1)?[thet3flag]:[[clotflag],[thet3flag]]
               doneclot[J]=1 
               if (verbose gt 2) then print,format='(%"found triplet %i %i %i")',sign*(*bundle[iClot,1])[0],sign*(*bundle[iClot,1])[1],sign*(*bundle[iClot,1])[2]
               break
            endif
         endfor
      endfor
      ;; if we cannot fill the complete configuration, drop it.
      if total(doneclot) ne n_elements(doneclot) then begin
         if (verbose gt 1) then print,format='(%"Data for time %f has only %d closures that include reference telecope %d, dropped.")',*bundle[0,0],total(doneclot),mytel1
         if (free_numberoftels) then warnHasUncompleteConfigurationSet++
         goto, FAILED_TO_ADD_STRUCT
      endif

      ;; create an ordered baseline list 12 13 14 23 24 34
      mybases=intarr(2,nbases)
      index=0
      for J=0L,ntels-1 do begin 
         mytel1=mytels[j] 
         for K=J+1,ntels-1 do begin 
            mytel2=mytels[K] 
            mybases[*,index]=[mytel1,mytel2] 
            index++ 
         endfor 
      endfor
      ;; create an ordered ucoord list
      mycoord=intarr(2,ncoords)
      index=0
      mytel1=mytels[0]
      for J=1L,ntels-1 do begin 
         mytel2=mytels[j] 
         mycoord[*,index]=[mytel1,mytel2] 
         index++ 
      endfor
      if (verbose gt 2) then print,'Closure telescope list used:'
      if (verbose gt 2) then print, myclos
      if (verbose gt 2) then print,'Baseline list used:'
      if (verbose gt 2) then print, mybases
      if (verbose gt 2) then print,'uv coordinates list used:'
      if (verbose gt 2) then print, mycoord
      ;;find and populate output struct wave
      wlen=*bundle[0,11]
      nbwlen=n_elements(wlen)
      donebase=intarr(nbases)
      for J=0L,nbases-1 do begin
         if (donebase[J] EQ 1) then break
         oucode1=oifits_baseline_encode(mybases[*,J],0,1)
         oucode2=oifits_baseline_encode(mybases[*,J],1,0)
         if (verbose gt 2) then print,format='(%"looking for base %d %d code %s or %s")',mybases[0,J],mybases[1,J],oucode1,oucode2
         for iClot=0L,nclot-1 do begin
            if (donebase[J] EQ 1) then break
            tripletbases=*bundle[iClot,10]
            for K=0,2 do begin
               icode=oifits_baseline_encode(*tripletbases[K],0,1)
               if ( icode EQ oucode1 OR icode EQ oucode2 ) then begin
                  if (verbose gt 2) then print,format='(%"found code %s at base %d of closure %d")',oucode1,K,iClot
                  thevis2=[*(*bundle[iClot,5])[K]]
                  thevis2err=[*(*bundle[iClot,6])[K]]
                  thevis2flag=*(*bundle[iClot,7])[K]
                  www=where(thevis2flag eq 70, count) & if (count GT 0) then thevis2flag[www]=0
                  www=where(thevis2flag eq 84, count) & if (count GT 0) then thevis2flag[www]=1

                  wherenan=where(~finite(thevis2), nancount) & if(nancount GT 0) then thevis2flag[wherenan]=1
                  wherenan=where(~finite(thevis2err), nancount) & if(nancount GT 0) then thevis2flag[wherenan]=1

                  vis2err=(total(size(vis2err)) LT 1)?[thevis2err]:[[vis2err],[thevis2err]]
                  vis2=(total(size(vis2)) LT 1)?[thevis2]:[[vis2],[thevis2]]
                  vis2flag=(total(size(vis2flag)) LT 1)?[thevis2flag]:[[vis2flag],[thevis2flag]]
                  if (hasOiVis) then begin
                     thevisphi=[*(*bundle[iClot,12])[K]]*!DPI/180.0
                     thevisphierr=[*(*bundle[iClot,13])[K]]*!DPI/180.0
                     thevisflag=*(*bundle[iClot,14])[K]
                     www=where(thevisflag eq 70, count) & if (count GT 0) then thevisflag[www]=0
                     www=where(thevisflag eq 84, count) & if (count GT 0) then thevisflag[www]=1

                     wherenan=where(~finite(thevisphi), nancount) & if(nancount GT 0) then thevisflag[wherenan]=1
                     wherenan=where(~finite(thevisphierr), nancount) & if(nancount GT 0) then thevisflag[wherenan]=1

                     visphi=(total(size(visphi)) LT 1)?[thevisphi]:[[visphi],[thevisphi]]
                     visphierr=(total(size(visphierr)) LT 1)?[thevisphierr]:[[visphierr],[thevisphierr]]
                     visflag=(total(size(visflag)) LT 1)?[thevisflag]:[[visflag],[thevisflag]]
                     if (hasOiVisData) then begin
                        thevisdata=[*(*bundle[iClot,15])[K]]
                        visdata=(total(size(visdata)) LT 1)?[thevisdata]:[[visdata],[thevisdata]]
                        theviserr=[*(*bundle[iClot,16])[K]]
                        viserr=(total(size(viserr)) LT 1)?[theviserr]:[[viserr],[theviserr]]
                     endif
                  endif
                  donebase[J]=1
                  break
               endif
            endfor 
         endfor
      endfor
      ;; vis2 order is 1stbase,nwlen,2nbase,nwlen,etc...
      donecoord=intarr(ncoords)
      for J=0L,ncoords-1 do begin
         if (donecoord[J] EQ 1) then break
         oucode1=oifits_baseline_encode(mycoord[*,J],0,1)
         oucode2=oifits_baseline_encode(mycoord[*,J],1,0)
         if (verbose gt 2) then print,format='(%"looking for base %d %d code %s or %s")',mycoord[0,J],mycoord[1,J],oucode1,oucode2
         for iClot=0L,nclot-1 do begin
            if (donecoord[J] EQ 1) then break
            tripletbases=*bundle[iClot,10]
            for K=0,2 do begin
               icode=oifits_baseline_encode(*tripletbases[K],0,1)
               if ( icode EQ  oucode1 ) then begin
                  if (verbose gt 2) then print,format='(%"found code %s at base %d of closure %d")',oucode1,K,iClot
                  thefrequ=(*bundle[iClot,8])[K]/*bundle[0,11]
                  freqs_u=(total(size(freqs_u)) LT 1)?[thefrequ]:[[freqs_u],[thefrequ]]
                  thefreqv=(*bundle[iClot,9])[K]/*bundle[0,11]
                  freqs_v=(total(size(freqs_v)) LT 1)?[thefreqv]:[[freqs_v],[thefreqv]]
                  donecoord[J]=1
                  break
               endif else if ( icode EQ  oucode2 ) then begin
                  if (verbose gt 2) then print,format='(%"found reverse code %s at base %d of closure %d")',oucode1,K,iClot
                  thefrequ=(*bundle[iClot,8])[K]/*bundle[0,11]
                  freqs_u=(total(size(freqs_u)) LT 1)?[-1.0D*thefrequ]:[[freqs_u],[-1.0D*thefrequ]]
                  thefreqv=(*bundle[iClot,9])[K]/*bundle[0,11]
                  freqs_v=(total(size(freqs_v)) LT 1)?[-1.0D*thefreqv]:[[freqs_v],[-1.0D*thefreqv]]
                  donecoord[J]=1
               endif 
            endfor 
         endfor
      endfor
      if nbWlen le 1 or keyword_set(noflatten) then begin
         vis2=transpose(vis2)&vis2err=transpose(vis2err)&vis2flag=transpose(vis2flag)&clot=transpose(clot)&cloterr=transpose(cloterr)
         clotflag=transpose(clotflag)&freqs_u=transpose(freqs_u)&freqs_v=transpose(freqs_v)&wlen=transpose(wlen) 
         struct={vis2:vis2*1D, vis2err:vis2err*1D, vis2flag:vis2flag, clot:clot*1D, cloterr:cloterr*1D, clotflag:clotflag, freqs_u:freqs_u*1D, freqs_v:freqs_v*1D, wlen:wlen*1D}
         if (hasOiVis) then begin
            visphi=transpose(visphi)&visphierr=transpose(visphierr)&visflag=transpose(visflag) 
            struct={vis2:vis2*1D, vis2err:vis2err*1D, vis2flag:vis2flag, clot:clot*1D, cloterr:cloterr*1D, clotflag:clotflag, freqs_u:freqs_u*1D, freqs_v:freqs_v*1D, wlen:wlen*1D, visphi:visphi*1D, visphierr:visphierr*1D, visflag:visflag }
         endif
         if (hasOiVisData) then begin
            visdata=transpose(visdata)&viserr=transpose(viserr)
            struct={vis2:vis2*1D, vis2err:vis2err*1D, vis2flag:vis2flag, clot:clot*1D, cloterr:cloterr*1D, clotflag:clotflag, freqs_u:freqs_u*1D, freqs_v:freqs_v*1D, wlen:wlen*1D, visphi:visphi*1D, visphierr:visphierr*1D, visdata:visdata, viserr:viserr, visflag:visflag}
         endif
         if (total(size(outdata)) LT 1) then begin
            desiredSize=ntimes
            outdata=replicate(struct,desiredSize)
         endif else begin
            if total(size(outdata[0].clot,/dim) ne  size(struct[0].clot,/dim)) gt 0 then begin
               if (free_numberoftels) then warnHasIncompatibleConfigurations++
               goto, FAILED_TO_ADD_STRUCT
            endif
            outdata[strindex]=struct
         endelse
         strindex++
      endif else begin
         vis2=transpose(reform(vis2, nbWlen, nbases,/overwrite))
         vis2err=transpose(reform(vis2err, nbWlen, nbases,/overwrite))
         vis2flag=transpose(reform(vis2flag, nbWlen, nbases,/overwrite))
         if (hasOiVis) then  visphi=transpose(reform(visphi, nbWlen, nbases,/overwrite))
         if (hasOiVis) then  visphierr=transpose(reform(visphierr, nbWlen, nbases,/overwrite))
         if (hasOiVis) then  visflag=transpose(reform(visflag, nbWlen, nbases,/overwrite))
         if (hasOiVisData) then visdata=transpose(reform(visdata, nbWlen, nbases,/overwrite))
         if (hasOiVisData) then viserr=transpose(reform(viserr, nbWlen, nbases,/overwrite))
         clot=transpose(reform(clot, nbWlen, nindeptclostheoric,/overwrite))
         cloterr=transpose(reform(cloterr, nbWlen, nindeptclostheoric,/overwrite))
         clotflag=transpose(reform(clotflag, nbWlen, nindeptclostheoric,/overwrite))
         freqs_u=transpose(reform(freqs_u, nbWlen, ncoords,/overwrite))
         freqs_v=transpose(reform(freqs_v, nbWlen, ncoords,/overwrite))
         if (hasOiVis) then begin
            for w=0L,nbWlen-1 do begin
               struct=(hasOiVisData)?{vis2:vis2[*,w]*1D, vis2err:vis2err[*,w]*1D, vis2flag:vis2flag[*,w], clot:clot[*,w]*1D, cloterr:cloterr[*,w]*1D, clotflag:clotflag[*,w], freqs_u:freqs_u[*,w]*1D, freqs_v:freqs_v[*,w]*1D, wlen:wlen[w]*1D, visphi:visphi[*,w]*1D, visphierr:visphierr[*,w]*1D, visdata:visdata[*,w], viserr:viserr[*,w], visflag:visflag[*,w]}:{vis2:vis2[*,w]*1D, vis2err:vis2err[*,w]*1D, vis2flag:vis2flag[*,w], clot:clot[*,w]*1D, cloterr:cloterr[*,w]*1D, clotflag:clotflag[*,w], freqs_u:freqs_u[*,w]*1D, freqs_v:freqs_v[*,w]*1D, wlen:wlen[w]*1D, visphi:visphi[*,w]*1D, visphierr:visphierr[*,w]*1D, visflag:visflag[*,w]}
               if (total(size(outdata)) LT 1) then begin
                  desiredSize=ntimes*nbWlen
                  outdata=replicate(struct,desiredSize)
               endif else begin
                  if size(outdata[0].clot,/dim) ne  size(struct[0].clot,/dim) then begin
                     if (free_numberoftels) then warnHasIncompatibleConfigurations++
                     goto, FAILED_TO_ADD_STRUCT
                  endif
                  outdata[strindex]=struct
               endelse
               strindex++
            endfor
         endif else begin
            for w=0L,nbWlen-1 do begin
               struct={vis2:vis2[*,w]*1D, vis2err:vis2err[*,w]*1D, vis2flag:vis2flag[*,w], clot:clot[*,w]*1D, cloterr:cloterr[*,w]*1D, clotflag:clotflag[*,w], freqs_u:freqs_u[*,w]*1D, freqs_v:freqs_v[*,w]*1D, wlen:wlen[w]*1D}
               if (total(size(outdata)) LT 1) then begin
                  desiredSize=ntimes*nbWlen
                  outdata=replicate(struct,desiredSize)
               endif else begin
                  if size(outdata[0].clot,/dim) ne  size(struct[0].clot,/dim) then begin
                     print,"FAILED_TO_ADD_STRUCT"
                     if (free_numberoftels) then warnHasIncompatibleConfigurations++
                     goto, FAILED_TO_ADD_STRUCT
                  endif
                  outdata[strindex]=struct
               endelse
               strindex++
            endfor
         endelse
      endelse
FAILED_TO_ADD_STRUCT:
      delvarx,vis2,vis2err,clot,cloterr,freqs_u,freqs_v,vis2flag,clotflag,wlen,struct
      if (hasOiVis) then delvarx,visphi, visphierr, visflag
      if (hasOiVisData) then delvarx,visdata,viserr
NOT_GOOD_ENOUGH:
   endfor
   ;; check big problem: outdata does not exist
   if (n_elements(outdata) lt 1) then  begin
      Message,/INFO,"SEVERE: Failure to produce requested data structure. Check OIFITS consistency and options of command. Eventually, refer to jmmc user support"
      return,-1
   endif
   ;;check final size:
   if (strindex ne desiredSize) then begin
      warnHasUncompleteConfigurationSet++
      outdata=(temporary(outdata))[0:strindex-1]
   endif
   if warnHasIncompatibleConfigurations gt 0 then  begin
      Message,/informational,"Data had incompatible configurations. Possible actions:"
      Message,/informational,"- Check the messages relative to configurations that appear when using the /verbose option."
      Message,/informational,"- Use the '/explode_triplets' option (but you lose important phase information)."
      Message,/informational,"- Alleviate synchronicity requirements using the 'synchro=xx' option where xx is in seconds."
   endif
   if warnHasUncompleteConfigurationSet gt 0 then begin
      Message,/informational,"Data had uncomplete configuration set, i.e., less closures than the number expected given the number of telescopes simultaneously used. The corresponding data has been dropped. Possible actions: "
      Message,/informational,"- use the '/explode_triplets' option (but you lose some information)."
   endif
   finalsize=n_elements(outdata.clot)
   if  (nbWlen le 1 or keyword_set(noflatten)) and (initialsize ne finalsize)  then message,/informational,"Step 2: Finally Retained "+string(finalsize,format='(I)')+" closures out of "+string(initialsize,format='(I)')+"." 
   if ~(nbWlen le 1 or keyword_set(noflatten)) and (initialsize*nbWlen ne finalsize) then message,/informational,"Step 2: Finally Retained "+string(finalsize/nbWlen,format='(%"%d")')+" non-redundant closures out of "+string(initialsize,format='(%"%d")')+" (redundant?) initially present in data."
   return,outdata
end 
