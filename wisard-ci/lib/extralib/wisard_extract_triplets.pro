;+
; NAME:
;     wisard_extract_triplets
; PURPOSE:
;     starting from all the closures of oit3 (see below), find all
;     data relevant to these closures and build a multi-dimensional
;     structure of pointers to be sorted by finction wisard_build_data(). This
;     structure contains all the exploitable triplets (since they have
;     been completed with the associated v2) but may contain complete
;     configurations (say, all the triplets for 6 tel) as well as less
;     complete configurations (if some closures were missing in the
;     data). As such, it must be filtered by another procedure to
;     feed into WISARD which admits only complete configurations at
;     the moment. 
;
; CALLING SEQUENCE:
;     out=wisard_extract_triplets(oiarray, oitarget, oiwavelength, oivis, oivis2, oit3, targetname=targetname, verbose=verbose, explode_triplets=explode_triplets, SYNCHRONICITY=synchronicity)
;
; INPUTS:
;     oiarray, oitarget, oiwavelength, oivis, oivis2, oit3
;     as described in read_oidata3.pro
;
; OPTIONAL INPUTS:
;     VERBOSE=[0..2]   be more.. eh.. verbose as number increases 
;     TARGETNAME : to select one target (mandatory if more than one in
;                  input data.
;     SYNCHRONICITY=value where value is in seconds the max admissible
;     difference between, e.g., V2 and corresponding T3 to consider
;     that they have been observed simultaneously.
;     EXPLODE: do not keep array configuration (all the triplets
;     pertaining to the same observational time) when creating the out
;     structure. This desperate move enables all the data to be fed in
;     WISARD whatever the disparity of the configurations in the
;     files(s) but at the expense of having much less phase information
;     (all the data is "observed" by a 3-telescope interferometer).
;
; OUTPUTS:
;      a self-containing unnamed structs of pointers, to be fed to
;      wisard_build_data().
;
; RESTRICTIONS:
;       None.
;
; PROCEDURE:
;
; EXAMPLE:
;
; PROCEDURES USED:
;
; MODIFICATION HISTORY:
;     v0.0 2013    G. Duvert
;-

function oifits_triplet_encode,array
; we sort the array to have a unique signature for a triangle
  on_error,2
  sarray=sort(array)
  return,FLOOR((array[sarray[0]]+1)*1000000L+(array[sarray[1]]+1)*1000L+(array[sarray[2]]+1)*1L)
end

function oifits_triplet_decode,value
; decode a triangle 
  on_error,2
  trois=floor(value/1000000L)
  deux=floor((value-trois*1000000L)/1000L)
  un=floor(value-trois*1000000L-deux*1000L)
  return,[trois-1,deux-1,un-1]
end
function oifits_triplet_compare,array1,array2
  on_error,2
; compare 2 triangles. Silly programmation. 
; same signature=+1 if same rotation -1 if inverted
; 0 otherwise (not equals)
  sign1=oifits_triplet_encode(array1)
  sign2=oifits_triplet_encode(array2)
  if sign1 ne sign2 then return,0
; triplets are identical -- check rotation
  plus=[[0,1,2],[1,2,0],[2,0,1]]
  s1=max(sort(array1)#plus) eq 5
  s2=max(sort(array2)#plus) eq 5
  if (s1 eq s2 ) then return,1 else return,-1
end

;function oifits_triplet_compare,array1,array2
;  on_error,2
;; compare 2 triangles. Silly programmation. 
;; same signature=+1 if same rotation -1 if inverted
;; 0 otherwise (not equals)
;  sign1=oifits_triplet_encode(array1)
;; build all possibilities of 3 values
;  possibilities=[[0,1,2],[1,2,0],[2,0,1],[0,2,1],[1,0,2],[2,1,0]]
;  ok=intarr(6)
;  for ipos=0,5 do begin  
;     sign2=oifits_triplet_encode(array2[possibilities[*,ipos]])
;     ok[ipos]= (sign1 eq sign2)
;  endfor
;  if total(ok) eq 0 then return,0
;; triplets are identical -- check rotation.
;  plus=[[0,1,2],[1,2,0],[2,0,1]]
;  s1=max(sort(array1)#plus) eq 5
;  s2=max(sort(array2)#plus) eq 5
;  if (s1 eq s2 ) then return,1 else return,-1
;end


function oifits_baseline_encode,array,i,j
  on_error,2
  return,FLOOR((array[i]+1)*1000L+(array[j]+1)*1L)
end

function oifits_baseline_decode,value
  on_error,2
  deux=floor(value/1000L)
  un=floor(value-deux*1000L)
  return,[deux-1,un-1]
end

function oifits_baseline_compare,array1,array2
 ON_ERROR, 2
 s11=oifits_baseline_encode(array1,0,1)
 s21=oifits_baseline_encode(array2,0,1)
 s22=oifits_baseline_encode(array2,1,0)
 if ((s11 eq s21) or (s11 eq s22)) then return,1 else return,0
end

function wisard_extract_triplets, oiarray, oitarget, oiwavelength, oivis, oivis2, oit3, targetname=targetname, verbose=verbose, explode_triplets=explode_triplets, SYNCHRONICITY=synchronicity
  ON_ERROR, 2

@ "wisard_common.pro"
@ "wisard_catch_noniteractive.pro" ; for interactive & catch facility.

; definitions of structures
; t3
t3time_=0
t3mjd_=1
t3int_time_=2
t3u1coord_=3
t3v1coord_=4
t3u2coord_=5
t3v2coord_=6
t3sta_index_=7
t3amp_=8
t3amperr_=9
t3phi_=10
t3phierr_=11
t3flag_=12
; v2 
v2time_      =0
v2mjd_       =1
v2int_time_  =2
v2ucoord_    =3
v2vcoord_    =4
v2sta_index_ =5
v2vis2data_  =6
v2vis2err_   =7
v2flag_      =8
;vis
vistime_     =0
vismjd_      =1
visint_time_ =2
visucoord_   =3
visvcoord_   =4
vissta_index_=5
visamp_   =6
visamperr_=7
visphi_   =8
visphierr_=9
visflag_=10
visdata_=11
viserr_ =12
;output structure
outtime_=0
outsta_index_=1
outt3phi_=2
outt3phierr_=3
outt3flag_=4
outv2_=5
outv2err_=6
outv2flag_=7
outucoord_=8
outvcoord_=9
outv2telindex_=10
outwl_=11
outvisphi_=12
outvisphierr_=13
outvisflag_=14
outvisdata_=15
outviserr_=16
;
secondInDays=1d/86400d

  if (n_elements(verbose) eq 0) then verbose=0
  if (n_elements(synchronicity) eq 0) then synchronicity=0.0d

; nb of objects:
  nbTargets=n_elements(oitarget)
  if (n_elements(targetname) eq 0 and nbTargets gt 1 ) then begin
     Print,"List Of targets in this dataset: ",oitarget.target
     MESSAGE,"Fatal -- More than 1 target in data, please specify using targetname= option"
  endif
  if  (n_elements(targetname) eq 0) then begin
     targetname=strtrim(oitarget.target,2)
     Print,"Defaulting to target: ",targetname
  endif
  whereismytarget=where(strtrim(oitarget.target,2) EQ targetname, count)
  if (count le 0) then MESSAGE,"Fatal -- Target "+targetname+" not found in data!"
  if (count gt 1) then MESSAGE,"Fatal -- Internal Error, Target multiply defined in data, please report!"
  itarget=whereismytarget

  if (size(oiarray,/dim) EQ 0) THEN MESSAGE,"Fatal -- no Array defined in data"
  if (size(oitarget,/dim) EQ 0) THEN MESSAGE,"Fatal -- no Targets defined in data"
  if (size(oiwavelength,/dim) EQ 0) THEN MESSAGE,"Fatal -- no Wavelengths defined in data"
  if (size(oit3,/dim) EQ 0) THEN MESSAGE,"Fatal -- no Closures present in data"
  if (size(oivis2,/dim) EQ 0) THEN MESSAGE,"Fatal -- no Visibilities present in data"
; oivis present?
  if (size(oivis,/dim) EQ size(oivis2,/dim) ) then hasOiVis=1 else hasOiVis=0 ; only this case can be used yet
  hasOiVisData=0
  if n_tags(oivis) ne 0 then if ( total(strpos(tag_names(oivis),"VISDATA")+1) gt 0 ) then hasOiVisData=(n_elements(oivis.visdata) gt 0)

  if (size(oivis,/dim) NE 0 AND NOT hasOiVis ) then begin
     hasOiVisData=0
     message,/informational,"WARNING--Data contains Insufficient number of Differential Visibilities, Will *NOT* use them"
  endif
; data for targets

; for itarget=0,nbTargets-1 do begin
  itarget=whereismytarget
  target=oitarget[itarget].target_id
  targnam=oitarget[itarget].target
  wheretargetv2=where(oivis2.target_id EQ target) 
  nv2=n_elements(oivis2[wheretargetv2])
  wheretargetT3=where(oit3.target_id EQ target)
  nt3=n_elements(oit3[wheretargett3])
  print, format='($, A12," has ",I0," V2 ",I0," T3 ")', targnam ,nv2, nt3
  if (hasOiVis) then begin
     wheretargetv=where(oivis.target_id EQ target)
     nv=n_elements(oivis[wheretargetv])
     print, format='(I0," V.")', nv 
  endif

;  if nt3 gt nv2/3 then print,"Data has (apparently) redundants closures"
;  if nt3 lt nv2/3 then print,"Warning -- Data has less closures than V2, WISARD will not be able to use all information"
;  
; find all the wavelength concerned by this target. we need T3, so the
; best is to start with T3 which is smaller to explore.
  aa=sort(oit3[wheretargett3].insname) & t3sortedinsnames=oit3[aa].insname
  list_of_backends=oit3[uniq(t3sortedinsnames)].insname
  nbackends=n_elements(list_of_backends)
  print,format='(%"%s %d %s")','found: ',nbackends,' spectral setup for this object'
; we have nbackends, so our observations are decomposed in nbackends
; how many t3 per backend?
  
  perBackendT3List=ptrarr(nbackends)
  for J=0L,nbackends-1 do begin  
     perBackendT3List[j]=ptr_new(where(oit3[wheretargett3].insname eq list_of_backends[J]) )  
  endfor

  perBackendV2List=ptrarr(nbackends)
  for J=0L,nbackends-1 do begin   
     perBackendV2List[j]=ptr_new(where(oivis2[wheretargetv2].insname eq list_of_backends[J]) )  
  endfor

  if (hasOiVis) then begin
     perBackendVisList=ptrarr(nbackends)
     for J=0L,nbackends-1 do begin
        perBackendVisList[j]=ptr_new(where(oivis[wheretargetv].insname eq list_of_backends[J]) )
     endfor
  endif

  perBackendWlList=ptrarr(nbackends)
  for K=0L,nbackends-1 do begin    
     for j=0L,n_elements(oiwavelength)-1 do begin   
        if (oiwavelength[J].insname EQ list_of_backends[K]) then begin   
           perBackendWlList[K]=ptr_new([*(oiwavelength[J]).eff_wave])   
           break   
        endif   
     endfor  
  endfor

; Now big loop on backends:
; how many telescopes involved in closures? for each time, select all
; the closures at this time. Important: since TIME and MJD may be
; incorrect (null) we choose MJD, but test if time gives more 'times' than
; TIME, on which case we choose time, but add date_obs as to convert
; to a simili-MJD (since we may have the same time in other,
; unrelated, observations).

  for iBackend=0L,nbackends-1 do begin
     if verbose gt 0 then print,format='(%"data from instrument number %d name %s")',iBackend,list_of_backends[iBackend]
     theT3=*(perBackendT3List[iBackend])             & if total(theT3 eq -1) eq 1 then CONTINUE
     theV2=*(perBackendV2List[iBackend])             & if total(theV2 eq -1) eq 1 then CONTINUE
     if (hasOiVis) then begin 
        theVis=*(perBackendVisList[iBackend]) 
        if total(theVis eq -1) eq 1 then CONTINUE
     endif

     t3_timelist=[oit3[theT3].mjd] 
     aa=sort(t3_timelist) & t3sortedtimes=t3_timelist[aa]
     t3times=t3sortedtimes[uniq(t3sortedtimes)]
     nt3mjd=n_elements(t3times) & if verbose gt 1 then print,"nt3mjd=",nt3mjd
     t3_timelist=[oit3[theT3].time] 
     aa=sort(t3_timelist) & t3sortedtimes=t3_timelist[aa]
     t3times=t3sortedtimes[uniq(t3sortedtimes)]
     nt3times=n_elements(t3times)  & if verbose gt 1 then print,"nt3times=",nt3times

     goodTimeIndex=1L                ; index of good time (MJD) in structures.
     if ( (nt3times gt nt3mjd) or (total([oit3[theT3].mjd]) le 1.0d) ) then begin
        message,'inconsistent times (mjd vs. time), using TIME+DATE_OBS',/informational
        goodTimeIndex=0              ; index of time
                                ; convert *all* times in mjd, by adding
                                ; time/86400.0d to date_obs (already
                                ; in MJD)
        oit3[theT3].time=(oit3[theT3].time/86400.0d)+oit3[theT3].date_obs
        oivis2[theV2].time=(oivis2[theV2].time/86400.0d)+oivis2[theV2].date_obs
        if (HasOiVis) then oivis[theVis].time=(oivis[theVis].time/86400.0d)+oivis[theVis].date_obs
     endif
     ;recompute with correct time = goodTimeIndex
     t3_timelist=(goodTimeIndex eq 1L)?[oit3[theT3].mjd]:[oit3[theT3].time] 
     aa=sort(t3_timelist) & t3sortedtimes=t3_timelist[aa]
     t3times=t3sortedtimes[uniq(t3sortedtimes)]
     nt3times=n_elements(t3times)  & if verbose gt 1 then print,"adopted nt3times=",nt3times

; we can construct a vector of pointers, one for each time, and a
; vector of [number of t3s involved at this time]
     t3_nt3=LONARR(nt3times)
     t3_struct=ptrarr(nt3times)
; now we fill that list with the subsets of t3 that correspond to each
; chunk of time:
     for iTime=0L,nt3times-1 do begin 
        aa=where(t3_timelist EQ t3times[iTime]) 
        t3_nt3[iTime]=n_elements(aa) 
        t3_struct[iTime]=ptr_new(aa) 
     endfor
;then, for each chunk of time, group all the t3 and v2
;data:

     myT3Values=ptrarr(nt3times,13) ;time,mjd,int_time,u1coo,v1coo,u2coo,v2coo,[t1,t2,t3],t3amp,t3amperr,t3phi,t3phierr,t3flag
     for iTime=0L,nt3times-1 do begin 
        subset=oit3[theT3[*(t3_struct[iTime])]] 
        myT3Values[iTime,t3time_]=ptr_new( [subset.time]) 
        myT3Values[iTime,t3mjd_]=ptr_new( [subset.mjd]) 
        myT3Values[iTime,t3int_time_]=ptr_new( [subset.int_time]) 
        myT3Values[iTime,t3u1coord_]=ptr_new( [subset.u1coord]) 
        myT3Values[iTime,t3v1coord_]=ptr_new( [subset.v1coord]) 
        myT3Values[iTime,t3u2coord_]=ptr_new( [subset.u2coord]) 
        myT3Values[iTime,t3v2coord_]=ptr_new( [subset.v2coord]) 
        myT3Values[iTime,t3sta_index_]=ptr_new( [subset.sta_index]) 
        myT3Values[iTime,t3amp_]=ptr_new( [subset.t3amp]) 
        myT3Values[iTime,t3amperr_]=ptr_new( [subset.t3amperr]) 
        myT3Values[iTime,t3phi_]=ptr_new([subset.t3phi]) 
        myT3Values[iTime,t3phierr_]=ptr_new([subset.t3phierr]) 
        myT3Values[iTime,t3flag_]=ptr_new([subset.flag]) 
     endfor


     myV2Values=ptrarr(9)       ;time,mjd,int_time,ucoo,vcoo,[t1,t2],vis2data,vis2err,flag
     myV2Values[v2time_]=ptr_new([oivis2[theV2].time])
     myV2Values[v2mjd_]=ptr_new([oivis2[theV2].mjd])
     myV2Values[v2int_time_]=ptr_new([oivis2[theV2].int_time])
     myV2Values[v2ucoord_]=ptr_new([oivis2[theV2].ucoord])
     myV2Values[v2vcoord_]=ptr_new([oivis2[theV2].vcoord])
     myV2Values[v2sta_index_]=ptr_new([oivis2[theV2].sta_index])
     myV2Values[v2vis2data_]=ptr_new([oivis2[theV2].vis2data])
     myV2Values[v2vis2err_]=ptr_new([oivis2[theV2].vis2err])
     myV2Values[v2flag_]=ptr_new([oivis2[theV2].flag])
;     nv2=n_elements(*myV2Values[v2time_])


     if (HasOiVis) then begin
        myVisValues=(hasOiVisData)?ptrarr(15):ptrarr(13) ;time,mjd,int_time,ucoo,vcoo,indexcode12,visamp,visamperr,visphi,visphierr,flag,[visdata,viserr]
        myVisValues[vistime_]=ptr_new([oivis[theVis].time])
        myVisValues[vismjd_]=ptr_new([oivis[theVis].mjd])
        myVisValues[visint_time_]=ptr_new([oivis[theVis].int_time])
        myVisValues[visucoord_]=ptr_new([oivis[theVis].ucoord])
        myVisValues[visvcoord_]=ptr_new([oivis[theVis].vcoord])
        myVisValues[vissta_index_]=ptr_new([oivis[theVis].sta_index])
        myVisValues[visamp_]=ptr_new([oivis[theVis].visamp])
        myVisValues[visamperr_]=ptr_new([oivis[theVis].visamperr])
        myVisValues[visphi_]=ptr_new([oivis[theVis].visphi])
        myVisValues[visphierr_]=ptr_new([oivis[theVis].visphierr])
        myVisValues[visflag_]=ptr_new([oivis[theVis].flag])
        if(hasOiVisData) then myVisValues[visdata_]=ptr_new([oivis[theVis].visdata])
        if(hasOiVisData) then myVisValues[viserr_]=ptr_new([oivis[theVis].viserr])
;        nv=n_elements(*myVisValues[vistime_])
     endif
; now, sort and write out the data structure
     dataSize=(HasOiVis)?15:12
     if (hasOiVisData) then dataSize=17
     for iTime=0L,nt3times-1 do begin ;different times
        actual_t3_dim=0L
        myT3bundle=ptrarr(t3_nt3[iTime],dataSize) ;t3time,[t1,t2,t3],clot,cloterr,clotflag,[vis2time will be same],*vis2[3],*vis2err[3],*vis2flag[3],ucoord[3],vcoord[3],tel1[3],tel2[3],*wlen + [visphi, visphierr, visflag] + [visdata, viserr]
        t3time=(*myT3Values[iTime,goodTimeIndex])[0]
        if verbose gt 1 then print,format='(%"time %d: %20.12f")',iTime,t3time

        t3timeMin=(*myT3Values[iTime,goodTimeIndex])[0]-(*myT3Values[t3int_time_]*secondInDays)/2.
        t3timeMax=(*myT3Values[iTime,goodTimeIndex])[0]+(*myT3Values[t3int_time_]*secondInDays)/2.

        for it3=0L,t3_nt3[iTime]-1 do begin ;same time, all concerned triplets
           t3_12= oifits_baseline_encode((*myT3Values[iTime,t3sta_index_])[*,it3],0,1)
           t3_23= oifits_baseline_encode((*myT3Values[iTime,t3sta_index_])[*,it3],1,2)
           t3_13= oifits_baseline_encode((*myT3Values[iTime,t3sta_index_])[*,it3],0,2)
           u1coord =(*myT3Values[iTime,t3u1coord_])[it3]
           v1coord =(*myT3Values[iTime,t3v1coord_])[it3]
           u2coord =(*myT3Values[iTime,t3u2coord_])[it3]
           v2coord =(*myT3Values[iTime,t3v2coord_])[it3]
           v2found=[-1,-1,-1]
           v2=ptrarr(3) & v2err=ptrarr(3) & v2flag=ptrarr(3)  &v2TelIndex=ptrarr(3) &v2base=intarr(3) &ucoord=dblarr(3) &vcoord=dblarr(3)
           ; list of v2 with time as t3time: v2 time available segment (+/- v2dit/2) should intersect the t3time segment (+/- t3dit/2)
           ; we should estimate the mutual coverage and take the greater one in case of doubt.

           ; first pass: normal OIFITS, times are OK OR Synchronicity has been correctly given
           candidatev2=where( abs(*myV2Values[goodTimeIndex]-t3time) le synchronicity*secondInDays, cv2count)
           ; they are not ok, try to use DIT
           if (cv2count lt 3) then begin
              ; we now try to use the DIT interval trick
              candidatev2=where( ( ( (*myV2Values[goodTimeIndex]) ge t3timeMin ) and ( (*myV2Values[goodTimeIndex]) le t3timeMax )  ), cv2count)
              if (cv2count lt 3) then begin
                 if verbose gt 1 then print,format='(%"available number of V2 (%d) prevents to fully characterize triplet no %d at time %20.12f")',cv2count,it3,t3time
                 break             ; no closure with less than 3!
              endif
           end

           for icount=0L,cv2count-1 do begin
              www=where(v2found eq -1, xxcount) & if (xxcount EQ 0) then break 
              iv2=candidatev2[icount]
              v2_12 = oifits_baseline_encode((*myV2Values[v2sta_index_])[*,iv2],0,1)
              v2_21 = oifits_baseline_encode((*myV2Values[v2sta_index_])[*,iv2],1,0)
              vindex=-1
              if (v2_12 EQ t3_12 OR v2_21 EQ t3_12) then vindex=0    ; would need to have a minus sign for other
              if (v2_12 EQ t3_23 OR v2_21 EQ t3_23) then vindex=1    ; observables, V2 is symmetrical wrt baseline sign
              if (v2_12 EQ t3_13 OR v2_21 EQ t3_13) then vindex=2 
              if (vindex GE 0) then begin
                 if ( v2found[vindex] eq -1) then begin
                    v2found[vindex]=0
                    v2[vindex]=(*myV2Values[v2vis2data_])[iv2]
                    v2err[vindex]=(*myV2Values[v2vis2err_])[iv2]
                    v2flag[vindex]=(*myV2Values[v2flag_])[iv2]
                    ucoord[vindex]=(*myV2Values[v2ucoord_])[iv2]
                    vcoord[vindex]=(*myV2Values[v2vcoord_])[iv2]
                    v2telindex[vindex]=ptr_new((*myV2Values[v2sta_index_])[*,iv2])
                    v2base[vindex]=(*myV2Values[v2sta_index_])[iv2]
; sanity check:
                    if vindex eq 0 and synchronicity eq 0.0d and abs(ucoord[0]-(*myT3Values[iTime,t3u1coord_])[it3]) gt 1D then message,/info,"inconsistent ucoord 0 "+string(ucoord[0],format='(G)')+string( (*myT3Values[iTime,t3u1coord_])[it3],format='(G)')+", please report!"
                    if vindex eq 0 and synchronicity eq 0.0d and abs(vcoord[0]-(*myT3Values[iTime,t3v1coord_])[it3])  gt 1D then message,/info,"inconsistent vcoord 0"+string(vcoord[0],format='(G)')+string( (*myT3Values[iTime,t3v1coord_])[it3],format='(G)')+", please report!"
                    if vindex eq 1 and synchronicity eq 0.0d and abs(ucoord[1]-(*myT3Values[iTime,t3u2coord_])[it3])  gt 1D then message,/info,"inconsistent ucoord 1"+string(ucoord[1],format='(G)')+string( (*myT3Values[iTime,t3u2coord_])[it3],format='(G)')+", please report!"
                    if vindex eq 1 and synchronicity eq 0.0d and abs(vcoord[1]-(*myT3Values[iTime,t3v2coord_])[it3])  gt 1D then message,/info,"inconsistent vcoord 1"+string(vcoord[1],format='(G)')+string( (*myT3Values[iTime,t3v2coord_])[it3],format='(G)')+", please report!"
                 endif
              endif
           endfor
           if (HasOiVis) then begin
              visfound=[-1,-1,-1]
              visphi=ptrarr(3) &visphierr=ptrarr(3) &visflag=ptrarr(3) &visdata=ptrarr(3) &viserr=ptrarr(3) 
              ; list of vis with time as t3time:
              candidatevis=where( *myVisValues[goodTimeIndex] eq t3time, cviscount)
              if (cviscount le 0) then break ;
              for icount=0L,cviscount-1 do begin
                 www=where(visfound eq -1, xxcount) & if (xxcount EQ 0) then break 
                 ivis=candidatevis[icount]
                 vis_12 = oifits_baseline_encode((*myVisValues[vissta_index_])[*,ivis],0,1)
                 vis_21 = oifits_baseline_encode((*myVisValues[vissta_index_])[*,ivis],1,0)
                 visindex=-1
                 sign=1
                 if (vis_12 EQ t3_12) then begin visIndex=0 & sign=1.0  & endif
                 if (vis_21 EQ t3_12) then begin visIndex=0 & sign=-1.0 & endif ; needs to have a minus sign for this observable
                 if (vis_12 EQ t3_23) then begin visIndex=1 & sign=1.0  & endif
                 if (vis_21 EQ t3_23) then begin visIndex=1 & sign=-1.0 & endif ; needs to have a minus sign for this observable
                 if (vis_12 EQ t3_13) then begin visIndex=2 & sign=1.0  & endif
                 if (vis_21 EQ t3_13) then begin visIndex=2 & sign=-1.0 & endif ; needs to have a minus sign for this observable
                 if (visIndex GE 0) then begin
                    if ( visfound[visIndex] eq -1) then begin
                       visfound[visIndex]=0
                       visphi[visIndex]=(*myVisValues[visphi_])[ivis] & *(visphi[visIndex])*=sign 
                       visphierr[visIndex]=(*myVisValues[visphierr_])[ivis]
                       visflag[visIndex]=(*myVisValues[visflag_])[ivis]
                       if (hasOiVisData) then begin
                          visdata[visIndex]=(*myVisValues[visdata_])[ivis]
                          viserr[visIndex]=(*myVisValues[viserr_])[ivis]
                          if (sign eq -1) then visdata[visIndex]=conj(*(visdata[visIndex]))
                       endif

                    endif
                 endif
              endfor
           endif

           www=where(v2found eq -1, v2foundcount)
           if (v2foundcount EQ 0) then begin
              if verbose gt 1 then print, format='(%"found complete triplet %d of %d (time %20.12f)")',it3+1,t3_nt3[iTime],t3time
              myT3bundle[actual_t3_dim,outtime_]=ptr_new(t3time)
              myT3bundle[actual_t3_dim,outsta_index_]=ptr_new((*myT3Values[iTime,t3sta_index_])[*,it3])            ; sta_index
              myT3bundle[actual_t3_dim,outt3phi_]=ptr_new((*myT3Values[iTime,t3phi_])[it3])             ; t3phi
              myT3bundle[actual_t3_dim,outt3phierr_]=ptr_new((*myT3Values[iTime,t3phierr_])[it3])             ; th3phierr
              myT3bundle[actual_t3_dim,outt3flag_]=ptr_new((*myT3Values[iTime,t3flag_])[it3])             ; t3flag
              myT3bundle[actual_t3_dim,outv2_]=ptr_new(v2)
              myT3bundle[actual_t3_dim,outv2err_]=ptr_new(v2err)
              myT3bundle[actual_t3_dim,outv2flag_]=ptr_new(v2flag)
              myT3bundle[actual_t3_dim,outucoord_]=ptr_new(ucoord)
              myT3bundle[actual_t3_dim,outvcoord_]=ptr_new(vcoord)
              myT3bundle[actual_t3_dim,outv2telindex_]=ptr_new(v2telindex)
              myT3bundle[actual_t3_dim,outwl_]=perBackendWlList[iBackend]

              if verbose gt 2 then begin
                 print, format='(%"t3time %20.12f\n")',t3time
                 print,(*myT3Values[iTime,t3sta_index_])[0:1,it3] , (*v2telindex[0]), ucoord[0], vcoord[0]
                 print,(*myT3Values[iTime,t3sta_index_])[1:2,it3] , (*v2telindex[1]), ucoord[1], vcoord[1]
                 print,(*myT3Values[iTime,t3sta_index_])[[0,2],it3] , (*v2telindex[2]), ucoord[2], vcoord[2]
              endif
; if found, add differential as well 
              if (HasOiVis) then begin
                 myT3bundle[actual_t3_dim,outvisphi_]=ptr_new(visphi)
                 myT3bundle[actual_t3_dim,outvisphierr_]=ptr_new(visphierr)
                 myT3bundle[actual_t3_dim,outvisflag_]=ptr_new(visflag)
                 if (hasOiVisData) then myT3bundle[actual_t3_dim,outvisdata_]=ptr_new(visdata)
                 if (hasOiVisData) then myT3bundle[actual_t3_dim,outviserr_]=ptr_new(viserr)
              endif
              actual_t3_dim++
           endif else begin
              if verbose gt 0 then print,format='(%"not enough V2 to fully characterize triplet no %d at time %20.12f")',it3,t3time
              goto, NO_T3_HERE
           endelse
        endfor
        if (actual_t3_dim  lt 1) then goto, NO_T3_HERE ; which can happen..
        ; create and populate output vis2,vis2err,vis2flag,clot,cloterr,clotflag,freqs_u,freqs_v
        ; do we have redundant closures (same time, same telescopes?)
        testarray=intarr(actual_t3_dim)
        for it3=0L,actual_t3_dim-1 do testarray[it3]=oifits_triplet_encode(*myT3bundle[it3,outsta_index_])
        test_actual_t3_dim=n_elements(uniq(testarray[sort(testarray)]))
        if test_actual_t3_dim NE actual_t3_dim then begin ; we have to reduce the number of triplets
           if verbose gt 0 then print, 'Warning, redundant triplets (same bases, same times) were present in data!'
           myT3bundle=myT3bundle[uniq(testarray[sort(testarray)]),*]
           actual_t3_dim=test_actual_t3_dim
        endif
        if (keyword_set(explode_triplets)) then begin
           if (total(size(outStructure)) eq 0 ) then begin
              outstructure=[ptr_new(myT3bundle[0,*])]
              for it3=1L,actual_t3_dim-1 do outstructure=[outstructure,ptr_new(myT3bundle[it3,*])]
           endif else begin
              for it3=0L,actual_t3_dim-1 do outstructure=[outstructure,ptr_new(myT3bundle[it3,*])]
           endelse
        endif else begin
           if (total(size(outStructure)) eq 0 ) then begin
              outstructure=[ptr_new(myT3bundle)] 
           endif else begin
              outstructure=[outstructure,ptr_new(myT3bundle)] 
           endelse
        endelse
NO_T3_HERE:
     endfor
  endfor
; endfor ; outer loop on targets supressed!
  if (total(size(outStructure)) eq 0) then begin 
     message,/informational,"Severe problem with this dataset, stopping here."
     message,/informational,"Not able to find at least one simultaneous complete (V2 and T3) observation in DataSet."
     message,"Consider using one of the options ""/merge_insnames"" or ""synchronicity=xx"" (xx in seconds) or edit the OI-FITS file. "
  endif
  size=0 & for i=0L,n_elements(outStructure)-1 do size+=(size(*outStructure[i],/dim))[0]
  if (size ne nt3) then message,/informational,"Step 1: Retained "+string(size,format='(I)')+" closures out of "+string(nt3,format='(I)')+" present initially in dataset."
  return, outStructure
end



