; $Id: wisard_jdata_x.pro,v 1.7 2008/01/08 16:42:03 mugnier Exp $ 
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
;   WISARD_JDATA_X - Generalised Mean square complex data likelihood
;                    Criterion for WISARD
;
; CATEGORY:
;	Array and Image Processing Routines 
;
; CALLING SEQUENCE:
;   value = WISARD_JDATA_X ( x, RE_PH = re_ph, IM_PH = im_ph
;                            , RE_Y_DATA = re_y_data, IM_Y_DATA = im_y_data 
;                            , W11 = w11, W12 = w12, W11 = w22
;                           [,GRADIENT_X=gradient_x] [, /VERSION] [, /HELP])
;
;
; PURPOSE:
;   This function computes and returns the value of the data likelihood
;   criterion of WISARD. Its gradient w.r.t. the object is
;   optionnally computed.
;
;
; POSITIONAL PARAMETERS:
;   x             : (input) 1D object for which we want to compute the value
;                   of the criterion.
;	
; KEYWORD PARAMETERS:
;   GRADIENT_X    : (optional input/output) if set and not equal to zero, contains in
;                    output the gradient w.r.t. to the object.
;   RE_PH         : (input) real part of matrix diag(exp i bbar alpha)#H
;   IM_PH         : (input) imaginary part of matrix diag(exp i bbar alpha)#H
;  
;   RE_Y_DATA     : (input) real part of data
;   IM_Y_DATA     : (input) imaginary part of data
;   INDX_FLAG     : (input) vectorial indices of "good"data
;  
;   W11,W12,W22   : (inputs) see doc scientifique
;
;   /VERSION : (optional input) prints version number before execution.
;   
;   /HELP    : (optional input) prints the documentation and exits.
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
; EXAMPLE:
;	Please provide a simple example here. 
;
; SEE ALSO:
;   All the WISARD_*.pro files, which are part of WISARD and covered by the
;   same copyright and license.
;   
; HISTORY:
;   $Log: wisard_jdata_x.pro,v $
;   
;   Revision 1.7  2008/01/08 16:42:03  mugnier
;   Fixed typo.
;
;   Revision 1.6  2007/10/31 12:48:36  mugnier
;   Added 3rd reference.
;
;   Revision 1.5  2007/01/29 10:38:29  meimon
;   *** empty log message ***
;
;   Revision 1.4  2006/10/31 18:27:38  meimon
;   cast� en double, grad_diff_fini OK
;
;   Revision 1.3  2006/10/27 18:28:21  meimon
;   Syntaxe renseign�e, compile
;
;   Revision 1.2  2006/10/27 17:37:08  meimon
;   *** empty log message ***
;
;   Revision 1.1  2006/10/27 17:13:03  meimon
;   Initial revision
;   
;-

FUNCTION  WISARD_JDATA_X, x, RE_PH = re_ph, IM_PH = im_ph, $
                          Ws=Ws, $
                          GRADIENT_X=gradient_x, $
                          VERSION = version, HELP = help, PRINT_TIMES=print_times

t_total=SYSTIME(/SECONDS )

;; To be changed manually if needed for optimization tests (=1)
print_times_loc=0

on_error,2
IF keyword_set(version) THEN $
    printf, -2, '% '+routine_courante()+': $Revision: 1.7 $, $Date: 2008/01/08 16:42:03 $'


IF (n_params() NE 1) OR keyword_set(help) THEN BEGIN ;(fill-in from syntax)
    message, 'Help required or incorrect syntax. Documentation:', /INFO
    doc_library, routine_courante()
    retall
ENDIF

if (print_times_loc) THEN t=SYSTIME(/SECONDS )

s=(size(re_ph))[1];
vect_w11 = (reform(*Ws.w11,s))
vect_w12 = (reform(*Ws.w12,s))
vect_w22 = (reform(*Ws.w22,s))

re_y_res = reform(*Ws.re_y_data,s)-re_ph#x
im_y_res = reform(*Ws.im_y_data,s)-im_ph#x

;if keyword_set(print_times) THEN print,'   >> wisard_jdata_x: re_ and im_y_res :',SYSTIME(/SECONDS )-t

;; computation of gradient, if defined and non-null : 
IF (gradient_x) THEN BEGIN  
 gradient_x_1 = matrix_multiply(vect_w11*re_y_res+vect_w12*im_y_res,re_ph,/ATRANSPOSE)
 gradient_x_2 = matrix_multiply(vect_w12*re_y_res+vect_w22*im_y_res,im_ph,/ATRANSPOSE) 
 gradient_x=-gradient_x_1-gradient_x_2
 ENDIF

IF (print_times_loc) THEN print,'   >> wisard_jdata_x: gradient_x :',SYSTIME(/SECONDS )-t
res=total(.5*vect_w11*abs2(re_y_res)+.5*vect_w22*abs2(im_y_res)+vect_w12*re_y_res*im_y_res)

if (print_times_loc) THEN print,'   >> wisard_jdata_x: TOTAL :',SYSTIME(/SECONDS )-t_total

return, res
END
