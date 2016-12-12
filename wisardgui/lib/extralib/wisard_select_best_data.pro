;+
; NAME:
;
; CATEGORY:
; Array and Image Processing Routines
;
; CALLING SEQUENCE:
; data = WISARD_SELECT_BEST_DATA( data )
;
; PURPOSE:
; When passed a 2-D data[i,j] structure, this program returns a 1d data
; structure where each structure of the first dimension of data (i) has
; been chosen according to some selection criterium between all the
; structures at the same [i] in the [j] possibilities
; Here, the criterium is the sum of inverses of w_tan and w_rad on cmdata,
; which is proportional to sigma^2.
;-



FUNCTION WISARD_SELECT_BEST_DATA, data, SIMULATED_DATA=SIMULATED_DATA
  on_error,2
  IF keyword_set(version) THEN $
     printf, -2, '% '+routine_courante()+': $Revision: 2.6 $, $Date: 2010-10-25 16:56:50 $'

  IF (n_params() NE 1) OR keyword_set(help) THEN BEGIN
     message, 'Help required or incorrect syntax. Documentation:', /INFO
     doc_library, routine_courante()
     retall
  ENDIF

  if ( (size(data))[0] eq 1 ) then return, data else begin
     goodData=data[*,0]                       ; define a good dataset, a priori.
     for ivariant=0,(size(data))[2]-1 do begin ; explore all datasets
        tmpdata=data[*,ivariant]
        n_tels = n_elements(tmpdata[0].freqs_u)+1
        operators = WISARD_OPERATORS(n_tels, VERSION = version)
        n_dat = n_elements(tmpdata)
        n_bases = operators.n_bases
        n_clot=(n_tels-1)*(n_tels-2)/2 ; not n_bases/3!
        tmpmdata = replicate({visamp:dblarr(n_bases), visamperr:dblarr(n_bases), $
                              visphi:dblarr(n_bases), visphierr:dblarr(n_bases), $
                              freqs_u:dblarr(n_bases), freqs_v:dblarr(n_bases)}, n_dat)
        tmpmdata.freqs_u = operators._B#tmpdata.freqs_u
        tmpmdata.freqs_v = operators._B#tmpdata.freqs_v
        tmpmdata.visamp=sqrt(0.5*(tmpdata.vis2+(tmpdata.vis2^2+2*tmpdata.vis2err^2)^0.5))
        tmpmdata.visamperr=1.0/sqrt(1.0/(tmpmdata.visamp^2)+2*(3*tmpmdata.visamp^2-tmpdata.vis2)/tmpdata.vis2err^2)
        tmpmdata.visphi = operators.dagC#tmpdata.clot ;;;Compute visphi
        if keyword_set(simulated_data) then begin
           tmpmdata.visphierr = 3.*((operators.dagC)^2#tmpdata.cloterr)
        endif else begin
           denom=operators.d#(0.5/(tmpmdata.visamp)^2)
           weightedCloterr2=tmpdata.cloterr^2/denom
           visphierr2 = 3.*((operators.dagC)^2#weightedCloterr2)/((tmpmdata.visamp)^2)
           tmpmdata.visphierr = sqrt(visphierr2)
        endelse
        ; make selection back at data level if noise is better
        if (ivariant eq 0) then cmdata=WISARD_MDATA2CMDATA(tmpmdata, OPERATORS = operators) else begin
           for j=0,n_elements(cmdata)-1 do begin
              tmpcmdata=WISARD_MDATA2CMDATA(tmpmdata, OPERATORS = operators)
              totold=total(sqrt(1.0/cmdata[j].w_rad^2+1.0/cmdata[j].w_tan^2))
              totnew=total(sqrt(1.0/tmpcmdata[j].w_rad^2+1.0/tmpcmdata[j].w_tan^2))
;;another possibility: favor long baselines closures.
;              totold=total(sqrt(cmdata[j].freqs_u^2+cmdata[j].freqs_v^2))
;              totnew=total(sqrt(tmpcmdata[j].freqs_u^2+tmpcmdata[j].freqs_v^2))
              if ( totnew lt totold ) then begin
;              if ( totnew gt totold ) then begin
                 print,format='(%"replacing %d (%f) -> (%f)")',j,totold,totnew
                 cmdata[j]=(WISARD_MDATA2CMDATA(tmpmdata,OPERATORS = operators))[j]
                 goodData[j]=tmpdata[j]
              end
           endfor
        endelse
     endfor
     return,  goodData
  endelse
END

