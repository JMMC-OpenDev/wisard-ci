; $Id: wisard_jdata_alpha.pro,v 2.0 2009/11/05 16:42:16 vannier & mugnier Exp $ 
; This file is part of the WISARD software.
;
; Copyright (c) S. Meimon and L. Mugnier, ONERA, 2004-2007.
;
; This software, WISARD, is made of wisard.pro and all the wisard_*.pro files
; under the same directory, as well as the example batch file wisard_batch.pro
; 
; This software is copyright ONERA.
; The authors are Serge Meimon and Laurent Mugnier. Latest version by M. Vannier
; E-mail: meimon at onera.fr, mugnier at onera.fr, mvannier@unice.fr
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
;   WISARD_JDATA_ALPHA - Generalised Mean square complex data likelihood
;                        Criterion for WISARD
;
; CATEGORY:
;	Array and Image Processing Routines 
;
; CALLING SEQUENCE:
;   critere = WISARD_JDATA_ALPHA(alpha, _B = _B, $
;                          ABS_HX = abs_hx, ARG_HX = arg_hx, $
;                          ABS_Y_DATA = abs_y_data, ARG_Y_DATA = arg_y_data, $
;                          WX1 = wx1, WX2 = wx2, W_RAD = w_rad, W_TAN = w_tan, $
;                          [, GRADIENT_ALPHA=gradient_alpha] $
;                          [, /VERSION] [, /HELP])
;
;
; PURPOSE:
;   This function computes and returns the value of the data likelihood
;   criterion of WISARD. Its gradient w.r.t. the aberrations alpha is
;   optionnally computed.
;
;
; POSITIONAL PARAMETERS:
;   alpha            : (input) aberrations 2D [n_tels-1 x n_nuplets] for which
;                      we compute the criterion
;	
; KEYWORD PARAMETERS:
;
;   GRADIENT_ALPHA   : (optional input/output) if set and not equal to zero, contains in
;                      output the gradient w.r.t. to the aberrations.
;   _B               : matrix Bbar
;   delta_alpha      : (output) Bbar#alpha
;   T_B              : transpos�e de la matrice Bbar
;
;   ABS_HX           : (input) array 2D [n_bases x n_nuplets] of |H#x|. 
;   ARG_HX           : (input) array 2D [n_bases x n_nuplets] of arg(H#x). 
;   ABS_Y_DATA       : (input) array 2D [n_bases x n_nuplets] of |y_data|. 
;   ARG_Y_DATA       : (input) array 2D [n_bases x n_nuplets] of arg(y_data). 
;
;   W_RAD            : (input) array 2D [n_bases x n_nuplets] of 1/var_rad. 
;   W_TAN            : (input) array 2D [n_bases x n_nuplets] of 1/var_tan. 
;
;   WX1              : (input) array 2D [n_bases x n_nuplets] of |H#x|^2(w_tan-w_rad). 
;   WX2              : (input) array 2D [n_bases x n_nuplets] of 2 w_rad|H#x| |y_data|. 
;
;
;   /VERSION : (optional input) prints version number before execution.
;   
;   /HELP    : (optional input) prints the documentation and exits.
;
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
;   $Log: wisard_jdata_alpha.pro,v $
;   
;   Revision 2.0  2009/11/09 vannier
;   Added the [indx_tag] indices to take into account only good-tagged data
;   
;   Revision 1.9  2008/01/08 16:42:16  mugnier
;   Fixed typo.
;
;   Revision 1.8  2007/10/31 12:49:05  mugnier
;   Added 3rd reference.
;
;   Revision 1.7  2007/01/29 10:37:25  meimon
;   *** empty log message ***
;
;   Revision 1.6  2007/01/26 19:35:07  meimon
;   doc added
;
;   Revision 1.5  2007/01/09 11:17:15  meimon
;   ??
;
;   Revision 1.4  2006/10/31 19:17:11  meimon
;   debug d'un facteur -0,5
;   refaire les calculs pour trouver le -
;   attention aux valeurs de H et B, qui impactent sur la pr�cision du calcul de
;   gradnum
;   cast� en double, OK par (my)grad_diff_finies sur BC O4'
;
;   Revision 1.3  2006/10/30 10:43:40  meimon
;   correc du call_procedure
;
;   Revision 1.2  2006/10/30 09:51:10  meimon
;   syntaxe renseign�e, compile ok
;
;   Revision 1.1  2006/10/27 18:15:49  meimon
;   Initial revision
;
;-




FUNCTION  WISARD_JDATA_ALPHA, alpha, DELTA_ALPHA = delta_alpha, _B = _B, T_B = t_B, $
                              ABS_HX = abs_hx, ARG_HX = arg_hx, $
                              ABS_Y_DATA = abs_y_data, ARG_Y_DATA = arg_y_data, $
                              WX1 = wx1, WX2 = wx2, W_RAD = w_rad, W_TAN = w_tan, $
                              GRADIENT_ALPHA = gradient_alpha, $
                              VERSION = version, HELP = help

on_error,2
IF keyword_set(version) THEN $
    printf, -2, '% '+routine_courante()+': $Revision: 1.9 $, $Date: 2008/01/08 16:42:16 $'

IF (n_params() NE 1) OR keyword_set(help) THEN BEGIN ;(fill-in from syntax)
    message, 'Help required or incorrect syntax. Documentation:', /INFO
    doc_library, routine_courante()
    retall
ENDIF


;IF NOT keyword_set(alpha) THEN message, 'alpha absents'

IF keyword_set(gradient_alpha) THEN gradient_alpha = alpha*0+1.

;calcul des r�sidus
delta_alpha = _B#alpha
phi_res = arg_y_data-arg_hx-delta_alpha
rad_y_res =  abs_y_data-abs_hx*COS(phi_res)
tan_y_res =  abs_hx*SIN(phi_res)
IF keyword_set(gradient_alpha) THEN gradient_alpha = -.5D*(t_B#(wx1*SIN(2*phi_res)+wx2*SIN(phi_res)))

return, total(.5D*w_rad*abs2(rad_y_res)+.5D*w_tan*abs2(tan_y_res))


END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;
;;%%%%%%%cr�ation des donn�es et de l'objet courant
;filename = '~/idl/BIC/PHASE1/final1.oifits'
;data = wisard_oifits2data(filename)
;unemas = 1d-3*(!DPi/180D)/3600D ; une mas en rd
;fov = 12.*unemas
;;seed = 1
;;ndat = 1000
;; data0 = {vis2:abs(randomn(seed, 15))*10., vis2err:abs(randomn(seed, 15)), clot:randomn(seed,10), $
;;          cloterr:abs(randomn(seed, 10)), freqs_u:abs(randomn(seed, 5))+1., freqs_v:abs(randomn(seed, 5))+1.}
;; data=replicate(data0,ndat)
;
;
;mdata = WISARD_DATA2MDATA(data, OPERATORS = operators)
;cmdata = WISARD_MDATA2CMDATA(mdata, OPERATORS = operators)
;H=wisard_make_h(FREQS_U=reform(cmdata.freqs_u, operators.n_bases*n_elements(mdata)),$ 
;                FREQS_V=reform(cmdata.freqs_v, operators.n_bases*n_elements(mdata)), $
;                FOV = fov, np_min = 64)
;H = H*1e10; pour le calcul de grad_num
;npix = (size(H))[2]
;sqnpix = sqrt(npix)
;TH = conj(transpose(H))
;dirty_map = TH#reform(cmdata.vis,(size(H))[1])
;xdm = randomn(seed, (size(H))[2])*10D;real_part(dirty_map) >0
;;aff, rotate(reform(xdm, sqnpix, sqnpix),5)
;alpha = randomn(seed, operators.n_tels-1, n_elements(cmdata))*1D
;
;
;re_y_data = real_part(cmdata.vis)
;im_y_data = imaginary(cmdata.vis)
;abs_y_data = abs(cmdata.vis)
;arg_y_data = atan(cmdata.vis,/phase)
;
;w11 = cmdata.w_rad*abs2(COS(arg_y_data)) $
;      +cmdata.w_tan*abs2(SIN(arg_y_data))
;w12 = cmdata.w_tan*abs2(COS(arg_y_data)) $
;      +cmdata.w_rad*abs2(SIN(arg_y_data))
;w22 = (cmdata.w_rad-cmdata.w_tan)*COS(arg_y_data)*SIN(arg_y_data)
;
;_B = operators._B; pour le calcul de grad_num
;t_B = transpose(_B)
;
;
;
;achix = reform(H#xdm, operators.n_bases, n_elements(cmdata))
;abs_hx = abs(achix)
;abs2_hx = abs2(achix)
;arg_hx = atan(achix,/phase)
;
;wx1 = abs2_hx*(cmdata.w_tan-cmdata.w_rad)
;wx2 = 2D*abs_hx*abs_y_data*cmdata.w_rad
;
;
;;%%%%%%%
;gradient = 1
;dummy = WISARD_JDATA_ALPHA(alpha, DELTA_ALPHA = delta_alpha, _B =_B, T_B = t_B, $
;                              ABS_HX = abs_hx, ARG_HX = arg_hx, $
;                              ABS_Y_DATA = abs_y_data, ARG_Y_DATA = arg_y_data, $
;                              WX1 = wx1, WX2 = wx2, W_RAD = cmdata.w_rad, W_TAN = cmdata.w_tan, $
;                           GRADIENT_ALPHA = gradient)
;
;P = EXP(!dI*reform(delta_alpha, operators.n_bases*n_elements(cmdata)))
;ph = diag(P)#H
;re_ph = real_part(H)
;im_ph = imaginary(H)
;dummix = WISARD_JDATA_X(xdm , RE_PH = re_ph, IM_PH = im_ph, RE_Y_DATA = $
;                        re_y_data, IM_Y_DATA = im_y_data, W11 = w11, W12 = w12, W22 = w22)
;print,  dummix/dummy
;END
;
;epsilon =.000001D
;
;gradnum=my_grad_diff_finies('wisard_jdata_alpha',$
;                         alpha, _B =_B, T_B = t_B, $
;                         ABS_HX = abs_hx, ARG_HX = arg_hx, $
;                         ABS_Y_DATA = abs_y_data, ARG_Y_DATA = arg_y_data, $
;                         WX1 = wx1, WX2 = wx2, W_RAD = cmdata.w_rad, W_TAN = cmdata.w_tan, $
;                         epsilon = epsilon)
;print, 'Corr�lation entre les 2 gradients :', correlate(gradient,gradnum)
;print, 'Rapport d''�chelle :', stddev(gradient)/stddev(gradnum)
;
;info, gradient-gradnum
;info, gradient
;
;END
