; $Id: j_prior_support.pro,v 1.2 2010-10-25 16:56:50 mvannier Exp $ ;
FUNCTION J_PRIOR_TOTVAR, x, GRADIENT_O = gradient, VERSION = version,  HELP = help

;+
; NAME:
;   J_PRIOR_SUPPORT - White quadratic regularization for soft support and more
;   
; CATEGORY:
;   Mathematics Routines
;
; CALLING SEQUENCE:
;   criterion = J_PRIOR_SUPPORT(X [,/VERSION] [, /HELP])
;
; PURPOSE:
;   The expression for the criterion is:
;   J_prior_totvar(X) = \sum_{p,q} [ grad_x(p,q) ]
;
;   
; POSITIONAL PARAMETERS:
;
;   X            : (input) object for which the regularization criterion is to be 
;                  computed. 
;   
;   MU           : (input) global hyper-parameter
;
;   /VERSION     : (input) prints the version number before execution.
;   
;   /HELP        : (input) prints the syntax and exits.
;
;
; SEE ALSO:
;   Other regularization criteria, notably  J_PRIOR_GAUSS and J_PRIOR_L1L2.
;
; AUTHOR :
; HISTORY:
;   Initial revision
;
;
;-

on_error,2
IF keyword_set(version) THEN $
    printf, -2, '% '+routine_courante()+': $Revision: 1.2 $, $Date: 2010-10-25 16:56:50 $'

IF (n_params() NE 1) OR keyword_set(help) THEN BEGIN ;(remplir d'apres syntaxe)
    message, 'Aide demandée ou syntaxe incorrecte. Documentation :', /INFO
    doc_library, routine_courante()
    retall
ENDIF
a=x-shift(x,0,1)
b=x-shift(x,0,-1)
c=x-shift(x,1,0)
d=x-shift(x,-1,0)
y=abs(a)+abs(b)+abs(c)+abs(d)
; gradient is sign1+sign2+sign3+sign4...
gradient=((a gt 0) - (a lt 0))+((b gt 0) - (b lt 0))+((c gt 0) - (c lt 0))+((d gt 0) - (d lt 0))
;wset,0
;tv,gradient
return, 0.25*total(y)
END

