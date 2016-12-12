function concat_oitable, tables, newtable

;+
; NAME:
;     concat_oitable
; PURPOSE:
;     In order to allow use of multiple OI-TABLES when using IDL
;     version 5.3 or less, one has to use a workaround to the 
;     standard concatenation methods.  This routine should allow the
;     OI_DATA IDL routines to be used on IDL version 5.1 or later 
;     (pointers must be supported!)
; 
;
; CALLING SEQUENCE:
;     combined_table = concat_oitable (tables,newtable)
;
; INPUTS:
;     tables : Existing array of structures
;   newtable : New table to add to 'tables'. Must be compatible.
;
; OPTIONAL INPUTS:
;     None 
;
; OUTPUTS:
;     Returns an array of structures, equivalent to [tables, newtable] under
;   	IDL 5.4.
;
; RESTRICTIONS:
;       None.
;
; PROCEDURE:
;      There are currently 6 different binary tables defined in the OIDATA
;      format.  This is a helper utility used to make arrays of the
;      structures used to represent the OI Tables.
;
;      The Data Exchange Standard for Optical (Visible/IR) Interferometry
;      is being maintained by the IAU Working group on Optical Interferometry
;      and the current version of the standard can be found at :
;            http://www.mrao.cam.ac.uk/~jsy1001/exchange/

; EXAMPLE:
;       This procedure is not meant to be called from the commandline.

; PROCEDURES USED:
;	This routine is part of the IDL Optical Interferometry (IOI) Library.
;        (more information at www.astro.lsa.umich.edu/~monnier)
;       The routines of the IOI Library generically 
;       require an up-to-date Astrolib library in order to read and write binary
;       FITS tables, among other things. The IDL Astronomy User's Library can
;       be found at http://idlastro.gsfc.nasa.gov/homepage.html.
;	
;
; MODIFICATION HISTORY:
;     v0.0 2003Feb18    J.D. Monnier    Written to extend functionality of
;	   library back to IDL 5.1
;     v0.1 2003Jull01	Takes advantage of newer IDL versions to speed up.
;-

IDL_VERSION = !version.RELEASE

if (idl_version le 5.1) then begin
a=tables[0]
c=tables
num=n_elements(newtable)

for i=0,num-1 do begin
 b=a
 struct_assign, newtable[i],b

 c=[c,b]
endfor
endif else begin
 c=[tables,newtable]
endelse

return,c
end

