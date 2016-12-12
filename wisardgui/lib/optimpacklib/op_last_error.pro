function op_last_error
;+
; NAME:
;   op_last_error
;
;
; PURPOSE:
;   Fetch last error message from OptimPack library.
;
;
; CATEGORY:
;   OptimPack.
;
;
; CALLING SEQUENCE:
;   op_last_error()
;
;
; INPUTS:
;   None.
;
;
; OPTIONAL INPUTS:
;   None.
;
;
; KEYWORD PARAMETERS:
;   None.
;
;
; OUTPUTS:
;   Last error message from OptimPack library is returned as a text string.
;   If there is no error message, the result is "no error".
;
; OPTIONAL OUTPUTS:
;   None.
;
;
; COMMON BLOCKS:
;   OP_COMMON - used to store the name of the shared OptimPack library.
;
;
; SIDE EFFECTS:
;   None.
;
;
; RESTRICTIONS:
;   Shared OptimPack library must be loaded (see OP_INIT).
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;   2003, Eric THIEBAUT.
;   $Id: op_last_error.pro,v 1.2 2010-11-30 21:26:11 duvert Exp $
;   $Log: not supported by cvs2svn $
;-
  common op_common, libname
  on_error, 2
  return, CALL_EXTERNAL(libname, 'op_idl_last_error', /S_VALUE)
end
