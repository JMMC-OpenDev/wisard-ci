FUNCTION WISARD_OIFITS2DATA, list_of_files, select_best_data=select_best_data, _EXTRA=ex
@ "wisard_common.pro"
@ "wisard_catch_noniteractive.pro" ; for interactive & catch facility.
read_oidata3, list_of_files, oiarray, oitarget,oiwavelength,oivis,oivis2, oit3, _EXTRA=ex
a=wisard_extract_triplets(oiarray, oitarget,oiwavelength,oivis,oivis2, oit3, _EXTRA=ex)
data=wisard_build_data(a, _EXTRA=ex)
return,data
end
