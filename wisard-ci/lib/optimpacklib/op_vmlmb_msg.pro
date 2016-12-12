; $Id: op_vmlmb_msg.pro,v 1.2 2010-11-30 21:26:11 duvert Exp $
function op_vmlmb_msg, ws
;+
; NAME:
;   OP_VMLMB_MSG
;
;
; PURPOSE:
;   Get message from workspace array used in OP_VMLMB.
;
;
; CATEGORY:
;   OptimPack.
;
;
; CALLING SEQUENCE:
;   op_vmlmb_msg(ws)
;
;
; INPUTS:
;   WS - workspace array as returned by op_vmlmb_setup and used by op_vmlmb.
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
;   A string message.
;
;
; OPTIONAL OUTPUTS:
;  None.
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
;   WS must be valid, shared OptimPack library must be loaded (see OP_INIT).
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
; NB: Revision 1.1 is the original version by Eric Thiebaut.
; $Log: not supported by cvs2svn $
; Revision 1.2  2007/01/24 11:04:34  mugnier
; Corrected 2 syntax errors in routine.
;
; Revision 1.1  2007/01/24 11:02:52  mugnier
; Initial revision
;
;
;-
  common op_common, libname
  on_error, 2
  return, CALL_EXTERNAL(libname, 'op_idl_vmlmb_msg', size(ws), ws, /S_VALUE)
end
