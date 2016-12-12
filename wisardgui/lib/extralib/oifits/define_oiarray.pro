pro define_oiarray,oiarray_unit

;+
; NAME:
;     define_oiarray
; PURPOSE:
;     Define the OIARRAY IDL Structure to hold data from the OI-DATA
;     FITS format OI_ARRAY binary table.  Currently returns structure
;     based on most current OI-DATA Revision.  Ability to specify
;     Revision number might be useful.
;
; CALLING SEQUENCE:
;     define_oiarray, oiarray_unit
;
; INPUTS:
;     None
;
; OPTIONAL INPUTS:
;     None 
;
; OUTPUTS:
;      OIARRAY_UNIT  Will  contain information on the interferometer array 
;		 geometry. This structure can be read by WRITE_OIDATA
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
;	This routine is part of the IDL Optical Interferometry (IOI) Library
;	 (more information at www.astro.lsa.umich.edu/~monnier)
;       The routines of the IOI Library generically 
;       require an up-to-date Astrolib library in order to read and write binary
;       FITS tables, among other things. The IDL Astronomy User's Library can
;       be found at http://idlastro.gsfc.nasa.gov/homepage.html.
;	
;
; MODIFICATION HISTORY:
;     v0.0 2002Jul02	J.D. Monnier	Written for OI-DATA Revision 
;     v0.1 2003Feb07	J.D. monnier	Updated to reflect 2002Nov26 update
;	   of the draft standard (Still officially Revision 0)
;     v1.0 2003Apr07	J.D. Monnier	Updated Revision to 1.0
;	
;
;     To do: A.  Add Revision Keyword (once multiple revisions exist)  
;-

oiarray_unit = { $
  extver : 0    ,$
  oi_revn: 1  	,$
  arrname: " "	,$
  frame:  " "  	,$
  arrayx: 0d0	,$
  arrayy: 0d0	,$
  arrayz: 0d0	,$
  tel_name: " " ,$
  sta_name: " "	,$
  sta_index   : 0	,$
  diameter: 1.0	,$
  staxyz  : dblarr(3)   $
 	}

end
