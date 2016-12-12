FUNCTION WISARD_OIFITS2DATA, list_of_files, select_best_data=select_best_data, _EXTRA=ex
read_oidata3, list_of_files, oiarray, oitarget,oiwavelength,oivis,oivis2, oit3, _EXTRA=ex
a=wisard_extract_triplets(oiarray, oitarget,oiwavelength,oivis,oivis2, oit3, _EXTRA=ex)
data=wisard_build_data(a, _EXTRA=ex)
if (keyword_set(select_best_data)) then begin
; find number of tels:
   ntels=(size(data[0].freqs_u,/dim))[0]+1
   for i=1,ntels-1 do begin
      tmp=wisard_build_data(a, number=ntels, first=i, _EXTRA=ex)
      data=[[data],[tmp]]
   end
   dataGood=wisard_select_best_data(data)
   return, dataGood
endif 

return,data
end
