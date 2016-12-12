function multiply_phasors, phasors,C
;
; Phasors is a matrix of elements(k,l) = exp( a(k,l) + j b(j,k)) (complex vis.) and C the closure to phase matrix.
; This function returns  the triple product of elements a(k,l) when C(k,l)=1 and conj(a(k,l)) when C(k,l)=-1
dims_phasors=size(phasors)
ndims_phasors=size(phasors,/n_dimensions)
dim2=(ndims_phasors eq 1)?1:dims_phasors[ndims_phasors]
dims_C=size(C)
result=complexArr(dims_C[1],dim2)+1 ; [10,120]

for l=0,dim2-1 do begin &$ ; 120
  for k=0,dims_C[1]-1 do begin &$; 10
    for i=0,dims_C[2]-1 do begin &$; 15
    if (C[k,i] EQ 1) then result[k,l]=result[k,l]*phasors[i,l] &$
    if (C[k,i] EQ -1) then result[k,l]=result[k,l]*conj(phasors[i,l]) &$
   endfor &$
   endfor &$
endfor
return, result
end
