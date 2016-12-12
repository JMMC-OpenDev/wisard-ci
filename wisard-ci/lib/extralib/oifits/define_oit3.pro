pro define_oit3,oit3_unit,nwave=nwave
;+
; NAME:
;     define_oit3
; PURPOSE:
;     Define the OIT3 IDL Structure to hold data from the OI-DATA
;     FITS format OI_T3 binary table.  Currently returns structure
;     based on most current OI-DATA Revision.  Ability to specify
;     Revision number might be useful.
;
; CALLING SEQUENCE:
;     define_oit3, oit3_unit, nwave=nwave
;
; INPUTS:
;     NWAVE : Number of Wavelengths for each data entry
;	      (normally defined in OI_WAVELENGTH). 
;	      Default: 1
;
; OPTIONAL INPUTS:
;     None 
;
; OUTPUTS:
;      OIT3_UNIT  Will contain measurements of the interferometer complex
;		 triple product. This structure can be read by WRITE_OIDATA
;		 to write standard OI-DATA FITS files.
;
; RESTRICTIONS:
;       None.
;
; PROCEDURE:
;      There are currently 6 different binary tables defined in the OIDATA
;      format.  This is one of the 6 procedures used to define these 
;      standard structures.
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
;     v0.0 2002Jul02    J.D. Monnier    Written for OI-DATA Revision
;     v0.1 2003Feb07    J.D. monnier    Updated to reflect 2002Nov26 update
;          of the draft standard (Still officially Revision 0)
;          Adopted use of pointers to variable-length arrays (NWAVE)
;     v1.0 2003Apr07    J.D. Monnier    Updated Revision to 1.0
;
;     To do: A.  Add Revision Keyword (once multiple revisions exist)  
;-


if (n_elements(nwave) eq 0) then nwave=1

oit3_unit = { $
  extver: 0 	,$
  nwave: fix(nwave)	,$
  oi_revn: 1	,$
  date_obs: " "	,$
  arrname: " "	,$
  insname: " "	,$
  target_id: 0	,$
  time: 0d0	,$
  mjd: 0d0	,$
  int_time: 0d0	,$
  t3amp: ptr_new(dindgen(nwave)), $
  t3amperr : ptr_new(dindgen(nwave)), $
  t3phi: ptr_new(dindgen(nwave)), $
  t3phierr: ptr_new(dindgen(nwave)), $
  u1coord : 0d0	,$
  v1coord : 0d0	,$
  u2coord : 0d0 ,$
  v2coord : 0d0 ,$
  sta_index : [1,2,3],$
  flag   : ptr_new(bytarr(nwave))	$   ; NO LOGICAL (1 bit) VARIABLES IN IDL
  }

end
