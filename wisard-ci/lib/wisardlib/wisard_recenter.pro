function wisard_recenter,x
; implies x is dim >=2
p=max(x,location,/NAN)
dim=size(x)
nx=dim[1]
ny=dim[2]
pos = ARRAY_INDICES(x, location)
return,shift(x,nx/2-pos[0],ny/2-pos[1])
end
