function cmdata_clot_maxfrequency, cmdata
nvis=(size(cmdata.vis,/dim))[0]
if (size(cmdata.vis))[0] eq 2 then nobs=(size(cmdata.vis,/dim))[1] else nobs = 0
ntels=round((1+sqrt(8.*nvis))/2.)
nclot=(ntels-1)*(ntels-2)/2

if (nobs eq 0) then begin
   temp=dblArr(3,nclot)
   index=0
   for J=0L,ntels-3 do begin 
      for K=J+1,ntels-2 do begin 
         temp[0,index] = cmdata.freqs_u[J]^2+cmdata.freqs_v[J]^2
         u0=cmdata.freqs_u[J]
         v0=cmdata.freqs_v[J]
         temp[1,index] = cmdata.freqs_u[K]^2+cmdata.freqs_v[K]^2
         temp[2,index] = (cmdata.freqs_u[K]-u0)^2+(cmdata.freqs_v[K]-v0)^2
         index++
      endfor 
   endfor
   result=sqrt(max(temp,dim=1))
   return, result
endif

temp=dblArr(3,nclot,nobs)
index=0
   for J=0,ntels-3 do begin 
      for K=J+1,ntels-2 do begin 
         temp[0,index,*] = cmdata.freqs_u[J]^2+cmdata.freqs_v[J]^2
         u0=cmdata.freqs_u[J]
         v0=cmdata.freqs_v[J]
         temp[1,index,*] = cmdata.freqs_u[K]^2+cmdata.freqs_v[K]^2
         temp[2,index,*] = (cmdata.freqs_u[K]-u0)^2+(cmdata.freqs_v[K]-v0)^2
         index++
      endfor 
   endfor
result=sqrt(max(temp,dim=1))
result=reform(result,nclot*nobs)
return, result
end
