pro wisard_average_data,masterDataArray,binsize
  nd=n_elements(masterDataArray)
  
  iloop=0
  for dindex=nd-1,0,-1 do begin
     iloop++
     data=(*masterDataArray[dindex])

     nbtimes=n_elements(data)
     nbwlen=n_elements(data[0].wlen)
     struct={vis2:data[0].vis2[*,0], vis2err:data[0].vis2err[*,0], vis2flag:data[0].vis2flag[*,0], clot:data[0].clot[*,0], cloterr:data[0].cloterr[*,0], clotflag:data[0].clotflag[*,0], freqs_u:data[0].freqs_u[*,0], freqs_v:data[0].freqs_v[*,0], wlen:data[0].wlen[0]}
     nout=nbWlen/binsize
     if (nout lt 1) then break
     desiredSize=nout*nbTimes
     outdata=replicate(struct,desiredSize)
     strindex=0
     for itime=0,nbtimes-1 do begin
        min=0
        max=binsize-1
        for i=0,nout-1 do begin
           vis2=data[itime].vis2[*,min:max]
           vis2err=data[itime].vis2err[*,min:max]
           vis2flag=data[itime].vis2flag[*,min:max]
           clot=data[itime].clot[*,min:max]
           cloterr=data[itime].cloterr[*,min:max]
           clotflag=data[itime].clotflag[*,min:max]
           freqs_u=data[itime].freqs_u[*,min:max]
           freqs_v=data[itime].freqs_v[*,min:max]
           wlen=data[itime].wlen[min:max]
           mvis2=mean(vis2,dim=2)
           mvis2err=mean(vis2err,dim=2)
           mvis2flag=mean(vis2flag,dim=2)
           mfreqs_u=mean(freqs_u,dim=2)
           mfreqs_v=mean(freqs_v,dim=2)
           wvis2=where(mvis2flag gt 0, count)
           if (count) then mvis2flag[wvis2]=1
           mvis2flag=byte(mvis2flag)
           mclot=mean(clot,dim=2)
           mcloterr=mean(cloterr,dim=2)
           mclotflag=mean(clotflag,dim=2)
           wclot=where(mclotflag gt 0, count)
           if (count) then mclotflag[wclot]=1
           mclotflag=byte(mclotflag)
           mwlen=mean(wlen)
           min+=binsize
           max+=binsize
           outdata[strindex++]={vis2:mvis2, vis2err:mvis2err, vis2flag:mvis2flag, clot:mclot, cloterr:mcloterr, clotflag:mclotflag, freqs_u:mfreqs_u, freqs_v:mfreqs_v, wlen:mwlen}
        endfor
     endfor
     masterDataArray[dindex]=ptr_new(outdata)
  endfor
end
