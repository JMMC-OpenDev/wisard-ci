; $Id: wisard_3t_find_phi.pro,v 1.7 2010-10-25 16:56:50 mvannier Exp $ 
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
;   edited by J�r�me Idier, ISTE, London, 2008.;

;+
; NAME:
;	WISARD_3T_FIND_PHI - computes the optimal starting phase data in the 3T case
;
; CATEGORY:
;	Array and Image Processing Routines 
;
; CALLING SEQUENCE:
; 
;	result = WISARD_3T_FIND_PHI(struct, operators [, /VERSION] [, /HELP])
;
;
; PURPOSE:
;   This function computes and returns the optimal phase to add to the pseudo
;   phase data CdagBeta (minimizing a criterion close to WISARD's),
;   given a guess object, for 3T data 
;
;
; POSITIONAL PARAMETERS:
;	struct      : (input) structure with H.guess, i.e. the model complex
;                         visibilities with the guess object
;   operators   : (input) structure, see function WISARD_OPERATORS 
;
;
;	
; KEYWORD PARAMETERS:
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
;
;   opt_phi = wisard_3t_find_phi(struct, operators)
;
; SEE ALSO:
;   WISARD_DATA2MDATA specifically
;
;   All the WISARD_*.pro files, which are part of WISARD and covered by the
;   same copyright and license.
;   
; HISTORY:
;   $Log: not supported by cvs2svn $
;   Revision 1.4  2008/01/08 11:49:53  mugnier
;    Fixed typo.
;
;   Revision 1.3  2007/10/31 12:39:33  mugnier
;   Added 3rd reference.
;
;   Revision 1.2  2007/10/30 10:22:24  meimon
;   *** empty log message ***
;
;   Revision 1.1  2007/10/01 14:19:21  meimon
;   Initial revision
;
;-



FUNCTION wisard_3t_find_phi,struct, operators, VERSION = version, HELP = help
on_error,2
IF keyword_set(version) THEN $
    printf, -2, '% '+routine_courante()+': $Revision: 1.7 $, $Date: 2010-10-25 16:56:50 $'

IF (n_params() NE 2) OR keyword_set(help) THEN BEGIN ;(fill-in from syntax)
    message, 'Help required or incorrect syntax. Documentation:', /INFO
    doc_library, routine_courante()
    retall
 ENDIF
print,  '3T find phi'
range = !Pi
decalx = 0
decaly = 0
nb = 50L*1D

opt_phi = replicate({phi:dblarr(2), crit:fltarr(nb^2)}, n_elements(struct))
;Note
;MATRIX_MULTIPLY can also be used in place of the ## operator. 
;For example, A ## B is equivalent to MATRIX_MULTIPLY(B, A), 
;and A ## TRANSPOSE(B) is equivalent to MATRIX_MULTIPLY(B, A, /ATRANSPOSE).
;i1 = findgen(NB)##replicate(1, NB)
i1=MATRIX_MULTIPLY(replicate(1, NB),findgen(NB))
;j = reform(transpose(findgen(NB)##replicate(1, NB)),nb^2)
j=reform(transpose(i1),nb^2)
i = reform(i1,nb^2)

; ii = ((i-((NB-1)/2))*(2*range/(NB-1)))
; jj = ((j-((NB-1)/2))*(2*range/(NB-1)))
ii = ((i)*(2*range/(NB-1)))
jj = ((j)*(2*range/(NB-1)))
ij = [[ii], [jj]]
tij = transpose(ij)
crit = fltarr(nb^2)
FOR i = 0, n_elements(struct)-1 DO BEGIN
   FOR j = 0, nb^2-1 DO BEGIN
      phi = tij[*, j]
      vec = (struct[i].resphi)- (operators._B#phi)
      ;crit[j] = transpose(vec)#diag((struct[i].resphierr)^(-2))#vec
      crit[j] = total(abs2(struct[i].visy*exp(!DI*operators._B#phi)-struct[i].achix))
   ENDFOR
;   aff, reform(crit, nb, nb)
   indmin = (where(crit EQ min(crit)))[0]
   opt_phi[i].crit = crit
   opt_phi[i].phi = tij[*, indmin] 
ENDFOR

return, opt_phi
END
