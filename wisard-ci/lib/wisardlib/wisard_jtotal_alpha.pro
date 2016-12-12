; $Id: wisard_jtotal_alpha.pro,v 1.6 2008/01/08 16:41:50 mugnier Exp $ 
; This file is part of the WISARD software.
;
; Copyright (c) S. Meimon and L. Mugnier, ONERA, 2004-2007.
;
; This software, WISARD, is made of wisard.pro and all the wisard_*.pro files
; under the same directory, as well as the example batch file wisard_batch.pro
; 
; This software is copyright ONERA.
; The authors are Serge Meimon and Laurent Mugnier.
; E-mail: meimon at onera.fr and mugnier at onera.fr
;
; This software is a computer program whose purpose is to reconstruct an image
; of an observed object from a set of interferometric measurements.
;
; This software is governed by the CeCILL-B license under French law and
; abiding by the rules of distribution of free software. You can use, modify
; and/ or redistribute the software under the terms of the CeCILL-B license as
; circulated by CEA, CNRS and INRIA at the following URL
; "http://www.cecill.info".
;
; As a counterpart to the access to the source code and  rights to copy,
; modify and redistribute granted by the license, users are provided only
; with a limited warranty  and the software's author,  the holder of the
; economic rights,  and the successive licensors  have only  limited
; liability. 
;
; In this respect, the user's attention is drawn to the risks associated
; with loading,  using,  modifying and/or developing or reproducing the
; software by the user in light of its specific status of free software,
; that may mean  that it is complicated to manipulate,  and  that  also
; therefore means  that it is reserved for developers  and  experienced
; professionals having in-depth computer knowledge. Users are therefore
; encouraged to load and test the software's suitability as regards their
; requirements in conditions enabling the security of their systems and/or 
; data to be ensured and,  more generally, to use and operate it in the 
; same conditions as regards security. 
;
; The fact that you are presently reading this means that you have had
; knowledge of the CeCILL-B license and that you accept its terms.
;
; Additionally, if you use this code or a code derived from it (e.g., for a
; publication), please cite the following papers, on which it is based:
;   
; (1) S. Meimon, L. M. Mugnier, and G. Le Besnerais,
;   "Reconstruction method for weak-phase optical interferometry", 
;   Opt. Lett., 30(14):1809-1811, July 2005.
;   PDF is on-line at: http://laurent.mugnier.free.fr/publis/Meimon-OL-05.pdf
; (2) S. Meimon, L. M. Mugnier, and Guy Le Besnerais,
;   "A convex approximation of the likelihood in optical interferometry",
;   J. Opt. Soc. Am. A, November 2005.
;   PDF is on-line at: http://laurent.mugnier.free.fr/publis/Meimon-JOSAA-05.pdf
; (3) L. M. Mugnier, G. Le Besnerais, and S. Meimon, 
;   "Inversion in optical imaging through atmospheric turbulence", 
;   chapter 10 of Bayesian Approach for Inverse Problems, 
;   edited by J�r�me Idier, ISTE, London, 2008.
;

;+
; NAME:
;	WISARD_JTOTAL_ALPHA - computes the criterion to minimize and its gradient
;                         w.r.t. the aberrations alpha
;
; CATEGORY:
;	Array and Image Processing Routines 
;
; CALLING SEQUENCE:
; value = WISARD_JTOTAL_ALPHA( alpha, gradient, _BALPHA =_Balpha
;                            ,WEIGHTS_CONSTANT = weights_constant
;                            , WEIGHTS_X = weights_x 
;                            , CRIT_ARRAY = crit_array 
;                            , OPERATORS = operators 
;                            [,/CHI2CRIT ][, /VERSION] [, /HELP])
;
;
;
; POSITIONAL PARAMETERS:
;	alpha:            (input ) current aberrations
;	gradient:         (optional input/output ) if set, the gradient is computed and
;                     stored in this keyword
;	
; KEYWORD PARAMETERS:
;
;  _BALPHA:           (input ) current differential baseline aberrations (see
;                     WISARD_OPERATORS.PRO)
;  WEIGHTS_CONSTANT:  (input ) see WISARD.PRO
;  WEIGHTS_X:         (input ) see WISARD.PRO
;  CRIT_ARRAY:        (input/output) contains the current value of the
;                                    criterion and the value of the
;                                    regularization terms 
;  OPERATORS:         (input ) see WISARD_OPERATORS.PRO
;
;  /CHI2CRIT : (optional input) casts the criterion in a chi2 style
;
;  /VERSION : (optional input) prints version number before execution.
;   
;  /HELP    : (optional input) prints the documentation and exits.
;
;
; AUTHORS:
; Serge Meimon and Laurent Mugnier (ONERA). 
; Contact: lastname at onera.fr
;
;
; RESTRICTIONS:
;   This code is copyright (c) S. Meimon and L. Mugnier, ONERA, 2004-2007. 
;   
;   WISARD is governed by the CeCILL-B license under French law and abiding by
;   the rules of distribution of free software. You can use, modify and/ or
;   redistribute the software under the terms of the CeCILL-B license as
;   circulated by CEA, CNRS and INRIA at the following URL:
;   "http://www.cecill.info".
;   See source code for the full notice.
;
;   Additionally, if you use this code or a code derived from it (e.g., for a
;   publication), please cite the following papers, on which it is based:
;   
;   (1) S. Meimon, L. M. Mugnier, and G. Le Besnerais,
;   "Reconstruction method for weak-phase optical interferometry", 
;   Opt. Lett., 30(14):1809-1811, July 2005.
;   PDF is on-line at: http://laurent.mugnier.free.fr/publis/Meimon-OL-05.pdf
;   (2) S. Meimon, L. M. Mugnier, and Guy Le Besnerais,
;   "A convex approximation of the likelihood in optical interferometry",
;   J. Opt. Soc. Am. A, November 2005.
;   PDF is on-line at: http://laurent.mugnier.free.fr/publis/Meimon-JOSAA-05.pdf
;   (3) L. M. Mugnier, G. Le Besnerais, and S. Meimon, 
;   "Inversion in optical imaging through atmospheric turbulence", 
;   chapter 10 of Bayesian Approach for Inverse Problems, 
;   edited by J�r�me Idier, ISTE, London, 2008.
;
;
; SEE ALSO:
;   All the WISARD_*.pro files, which are part of WISARD and covered by the
;   same copyright and license.
;   
; HISTORY:
;   $Log: wisard_jtotal_alpha.pro,v $
;   Revision 1.6  2008/01/08 16:41:50  mugnier
;   Fixed typo.
;
;   Revision 1.5  2007/10/31 12:48:11  mugnier
;   Added 3rd reference.
;
;   Revision 1.4  2007/10/30 10:20:52  meimon
;   chi^2 style crterion keyword added
;
;   Revision 1.3  2007/01/29 10:39:25  meimon
;   *** empty log message ***
;
;   Revision 1.2  2007/01/26 19:25:49  meimon
;   keyword method removed
;   doc added
;
;-

FUNCTION WISARD_JTOTAL_ALPHA, alpha, gradient, _BALPHA = _Balpha,WEIGHTS_CONSTANT = weights_constant, WEIGHTS_X = weights_x $
                              , CRIT_ARRAY = crit_array $
                              , OPERATORS = operators, CHI2CRIT = chi2crit, VERSION = version, HELP $
                              = help


;; put the following 9 lines (down to 1st ENDIF) just after routine name:
on_error,2
IF keyword_set(version) THEN $
    printf, -2, '% '+routine_courante()+': $Revision: 1.6 $, $Date: 2008/01/08 16:41:50 $'

IF (n_params() GT 2) OR keyword_set(help) THEN BEGIN ;(fill-in from syntax)
    message, 'Help required or incorrect syntax. Documentation:', /INFO
    doc_library, routine_courante()
    retall
ENDIF



IF (n_params() EQ 2) THEN gradient = 1. ELSE gradient = 0. 
crit_array[0] = WISARD_JDATA_ALPHA(alpha, _B = operators._B, DELTA_ALPHA = _Balpha, T_B = operators.t_B, $
                             ABS_HX = weights_x.abs_hx, ARG_HX = weights_x.arg_hx, $
                             WX1 = weights_x.wx1, WX2 = weights_x.wx2, $
                             ABS_Y_DATA = weights_constant.abs_y_data, ARG_Y_DATA = weights_constant.arg_y_data, $
                             W_RAD = weights_constant.w_rad, W_TAN = weights_constant.w_tan, $
                             GRADIENT_ALPHA = gradient, $
                             VERSION = version, HELP = help)
;print, (crit_array)
IF keyword_set(chi2crit) THEN mult = 1./n_elements(*weights_constant.re_y_data) ELSE mult = 1.
crit_array = crit_array*mult
gradient =  gradient*mult

;print, "wisard_jtotal_alpha total(crit_array): ";
;print, total(crit_array)
return, total(crit_array)

END
