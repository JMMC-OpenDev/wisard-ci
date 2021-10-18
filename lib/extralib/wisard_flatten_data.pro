pro wisard_flatten_data,masterDataArray
  nd=n_elements(masterDataArray)
  
  iloop=0
  for dindex=nd-1,0,-1 do begin
     iloop++
     data=(*masterDataArray[dindex])
     isnotflat= (size(data[0].wlen))[0] gt 0
     if (isnotflat) then begin
        nbtimes=n_elements(data)
        nbwlen=n_elements(data[0].wlen)
        struct={vis2:data[0].vis2[*,0]*1D, vis2err:data[0].vis2err[*,0]*1D, vis2flag:data[0].vis2flag[*,0], clot:data[0].clot[*,0]*1D, cloterr:data[0].cloterr[*,0]*1D, clotflag:data[0].clotflag[*,0], freqs_u:data[0].freqs_u[*,0]*1D, freqs_v:data[0].freqs_v[*,0]*1D, wlen:data[0].wlen[0]*1D}
        desiredSize=nbWlen*nbTimes
        outdata=replicate(struct,desiredSize)
        strindex=1
        for itime=0,nbtimes-1 do begin
           for w=1L,nbWlen-1 do begin
              outdata[strindex++]={vis2:data[itime].vis2[*,w]*1D, vis2err:data[itime].vis2err[*,w]*1D, vis2flag:data[itime].vis2flag[*,w], clot:data[itime].clot[*,w]*1D, cloterr:data[itime].cloterr[*,w]*1D, clotflag:data[itime].clotflag[*,w], freqs_u:data[itime].freqs_u[*,w]*1D, freqs_v:data[itime].freqs_v[*,w]*1D, wlen:data[itime].wlen[w]*1D}
           endfor
        endfor
        masterDataArray[dindex]=ptr_new(outdata)
     endif
  end
end
