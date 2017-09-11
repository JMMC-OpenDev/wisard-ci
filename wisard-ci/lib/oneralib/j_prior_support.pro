; $Id: j_prior_support.pro,v 1.2 2010-10-25 16:56:50 mvannier Exp $ ;
FUNCTION J_PRIOR_SUPPORT, x, GRADIENT_O = gradient, $
                          MU = mu, MEAN_O = mean_o, FWHM = fwhm, $
                          VERSION = version,  HELP = help

;+
; NAME:
;   J_PRIOR_SUPPORT - White quadratic regularization for soft support and more
;   
; CATEGORY:
;   Mathematics Routines
;
; CALLING SEQUENCE:
;   criterion = J_PRIOR_SUPPORT(X, GRADIENT_O = gradient, $
;                               MU = mu, $
;                               [MEAN_O = mean_o | FWHM = fwhm] $
;                               [, /VERSION] [, /HELP])
;
; PURPOSE:
;   J_PRIOR_SUPPORT implements a white quadratic regularization that draws the
;   reconstructed object X towards a given prior object MEAN_O.
;   This regularization was proposed by L. Mugnier (Onera) and E. Thiébaut
;   (CRAL) for the original purpose of implementing a soft support constraint, 
;   which is achieved by taking a Lorentzian as the prior object.
;   This regularization is specifically designed for Optical Interferometry
;   as MEAN_O is the mimnimizer of the criterion under the constraint of
;   unit sum of X.
;
;   The expression for the criterion is:
;   J_prior(X) = MU * \sum_{p,q} (X^2(p,q) / MEAN_O(p,q)
;
;   Reference:
;   Imaging with long-baseline optical interferometry, 
;   G. Le Besnerais, S. Lacour, L. M. Mugnier, E. Thiebaut, G. Perrin
;   and S. Meimon, IEEE JOURNAL OF SELECTED TOPICS IN SIGNAL PROCESSING,
;   submitted (2008). 
;   
; POSITIONAL PARAMETERS:
;
;   X            : (input) object for which the regularization criterion is to be 
;                  computed. 
;   
; KEYWORD PARAMETERS:
;
;   GRADIENT_O   : (input/output) if set to a non-zero variable, receives the
;                  gradient of the criterion on exit.
;
;   MU           : (input) global hyper-parameter
;
;   MEAN_O       : (optional input) prior object. If not set, the prior object
;                  is a Lorentzian defined by FWHM.
;
;   FWHM         : (optional input) FWHM, in pixels, of the Lorentzian used as
;                  default prior object if MEAN_O is not set. NB= this is
;                  costly because mean_o is re-computed at each call of this
;                  routine. 
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
;   $Author: mvannier $
;
; HISTORY:
;   From the routine j_support 1.2 by Serge Meimon (a priori = Lorentzian). 
;   $Log: not supported by cvs2svn $
;   Revision 1.2  2008/03/28 15:14:33  mugnier
;   Documentation and homogenization of keyword names.
;
;   Revision 1.1  2008/02/20 18:46:19  mugnier
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

IF NOT (keyword_set(mean_o) OR keyword_set(fwhm) ) THEN $
    message, 'il faut passer mean_o (prioritaire) ou fwhm'
np = sqrt(n_elements(x))

IF keyword_set(mean_o) THEN $
    prior_inside = mean_o $
ELSE BEGIN
   dist = dist(np, cx = np/2D, cy = np/2D)*1D;*1./np
   prior_inside = 1D / (1D + (2D*dist/fwhm)^2) ; Lorentzian shape by default
ENDELSE

; weights = mu / prior_inside

IF keyword_set(gradient) THEN gradient = (2D*mu) * (x/prior_inside)

return, mu * total(abs2(x)/prior_inside)

END


Np = 32
object= rebin(readfits('/home/opc/images/galaxie/m51-64x64.fits'), NP, NP, /SAM)

gradient = 1 & dummy = j_prior_support(object,gradient=gradient, fwhm = 3., $
                                       mu = 1D)
gradnum=grad_diff_finies('j_prior_support', object, fwhm = 3., mu = 1D, epsilon = .05);
; NB : il faut ajuster le epsilon a la main
print, 'Corrélation entre les 2 gradients :', correlate(gradient,gradnum)
print, 'Rapport d''échelle :', stddev(gradient)/stddev(gradnum)
aff, [gradient, gradnum, 1e2*(gradient - gradnum)]

END
