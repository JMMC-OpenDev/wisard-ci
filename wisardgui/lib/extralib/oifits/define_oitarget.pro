pro define_oitarget,oitarget_unit


;+
; NAME:
;     define_oitarget
; PURPOSE:
;     Define the OITARGET IDL Structure to hold data from the OI-DATA
;     FITS format OI_TARGET binary table.  Currently returns structure
;     based on most current OI-DATA Revision (see revision history at the
;     end of header).  Ability to specify Revision number might be useful.
;
; CALLING SEQUENCE:
;     define_oitarget, oitarget_unit
;
; INPUTS:
;     None
;
; OPTIONAL INPUTS:
;     None 
;
; OUTPUTS:
;      OITARGET_UNIT  Will  contain information on the interferometer 
;		 targets. This structure can be read by WRITE_OIDATA
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
;          of the draft standard (Still officially Revision 0). Note in
;	   this revision, there can not be multiple target tables (no 
;          EXTVER keyword needed).
;     v0.2 2003Feb13	J.D. Monnier	Updated to reflect modified
;	   column headings 
;     v1.0 2003Apr07    J.D. Monnier    Updated Revision to 1.0
;
;     To do: A.  Add Revision Keyword (once multiple revisions exist)  
;-

oitarget_unit = { $
  oi_revn : 1  	,$
  target_id:0	,$
  target: " "	,$
  raep0: 0d0	,$
  decep0: 0d0	,$
  equinox: 1.0	,$
  ra_err: 0d0	,$
  dec_err: 0d0	,$
  sysvel: 0d0	,$
  veltyp: " "	,$
  veldef: " "	,$
  pmra:	0d0	,$
  pmdec: 0d0    ,$
  pmra_err: 0d0	,$
  pmdec_err: 0d0,$
  parallax: 1.0 ,$
  para_err: 1.0	,$
  spectyp: " " 	$
     }

end
 
  
