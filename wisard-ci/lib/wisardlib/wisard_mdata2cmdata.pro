; $Id: wisard_mdata2cmdata.pro,v 1.11 2010-10-25 16:56:50 mvannier Exp $ 
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
;	WISARD_MDATA2CMDATA - Computes convexified myopic data from myopic data
;
; CATEGORY:
;	Array and Image Processing Routines 
;
; CALLING SEQUENCE:
; 
;	cmdata = WISARD_MDATA2CMDATA( mdata, OPERATORS = operators[, /DL][, /VERSION] [, /HELP])
;
;
; PURPOSE:
;   This function computes and returns convexified myopic data from myopic
;   data. The convexified myopic data are aranged in a array of Structures
;   containing the following fields (see doc scientifique):
;
;    VIS             DCOMPLEX  Array[N_BASES]
;    W_RAD           DOUBLE    Array[N_BASES]
;    W_TAN           DOUBLE    Array[N_BASES]
;    FREQS_U         DOUBLE    Array[N_BASES]
;    FREQS_V         DOUBLE    Array[N_BASES]
;
; POSITIONAL PARAMETERS:
;	mdata: (input) Array of Structures containing the following fields (see doc scientifique): 
;	 VISAMP          DOUBLE    Array[N_BASES]
;    VISAMPERR       DOUBLE    Array[N_BASES]
;    VISPHI          DOUBLE    Array[N_BASES]
;    VISPHIERR       DOUBLE    Array[N_BASES]
;    FREQS_U         DOUBLE    Array[N_BASES]
;    FREQS_V         DOUBLE    Array[N_BASES]
;
; KEYWORD PARAMETERS:
;
;   OPERATORS: (input) See WISARD_OPERATORS.
;
;   /DL      : (optional input) Uses the limited expansion of w_rad,w_tan.
;              Obsolete now. 
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
;	mdata = WISARD_DATA2MDATA(data, OPERATORS = operators, VERSION = version)
;  cmdata = WISARD_MDATA2CMDATA(mdata, OPERATORS = operators, VERSION = version)
;
; SEE ALSO:
;   All the WISARD_*.pro files, which are part of WISARD and covered by the
;   same copyright and license.
;   
; HISTORY:
;   $Log: not supported by cvs2svn $
;   Revision 1.9  2008/01/08 16:41:07  mugnier
;   Fixed typo.
;
;   Revision 1.8  2007/10/31 12:45:34  mugnier
;   Added 3rd reference.
;
;   Revision 1.7  2007/01/29 13:00:33  meimon
;   doc added
;
;   Revision 1.6  2007/01/25 14:32:47  meimon
;   bug sur fact_carre corrig�
;
;   Revision 1.5  2007/01/25 14:07:39  meimon
;   correction erreur de calcul de la these
;
;   Revision 1.4  2006/11/09 18:00:25  meimon
;   ajout keyword dl
;
;   Revision 1.3  2006/11/06 14:59:59  meimon
;   !I->!dI
;
;   Revision 1.2  2006/10/31 18:30:01  meimon
;   cast� en double, test sur BIC 04' ok
;
;   Revision 1.1  2006/10/30 14:37:22  meimon
;   Initial revision
;
;-

FUNCTION WISARD_MDATA2CMDATA, mdata, OPERATORS = operators, DL = dl, VERSION = version, HELP = help

on_error,2
IF keyword_set(version) THEN $
    printf, -2, '% '+routine_courante()+': $Revision: 1.11 $, $Date: 2010-10-25 16:56:50 $'

IF (n_params() NE 1) OR keyword_set(help) THEN BEGIN ;(fill-in from syntax)
    message, 'Help required or incorrect syntax. Documentation:', /INFO
    doc_library, routine_courante()
    retall
ENDIF


n_dat = n_elements(mdata)
IF NOT keyword_set(operators) THEN message, 'champ operators vide'
n_tels = operators.n_tels
n_bases = operators.n_bases


cmdata = replicate({vis:dcomplexarr(n_bases), w_rad:dblarr(n_bases),$
w_tan:dblarr(n_bases), freqs_u:dblarr(n_bases), freqs_v:dblarr(n_bases)}, $
    n_dat)

;;;Compute freqs
cmdata.freqs_u = mdata.freqs_u
cmdata.freqs_v = mdata.freqs_v


;;special test case with complex (radio, aspro2...) values
;if (n_tags(mdata) eq 4) then begin
;cmdata.vis = mdata.visdata
;visamp=sqrt(mdata.visdata^2)
;visphi=atan(mdata.visdata,/phase)
;visamperr=sqrt(mdata.viserr^2)
;visphierr=atan(mdata.viserr,/phase)
;fact_plus = .5D*(EXP(-2.D*visphierr^2)+1D)
;fact_moins = .5D*(-EXP(-2.D*visphierr^2)+1D)
;fact_carre = .5D*(-EXP(-visphierr^2)+1D)^2
;cmdata.w_rad = (fact_plus*abs2(visamperr)+fact_carre*abs2(visamp))^(-1)
;cmdata.w_tan = (fact_moins*(abs2(visamperr)+abs2(visamp)))^(-1)
;return, cmdata
;endif

;;;  Compute vis
; note GD: here we go away from "measured" if visphierr is large
cmdata.vis = mdata.visamp*EXP(!dI*mdata.visphi)*(-EXP(-.5D*mdata.visphierr^2)+2D)

;;;;; MY_TEST:
;gorb=(-EXP(-.5D*mdata.visphierr^2)+2D)
;print, 'facteur 2-exp : mean=',mean(gorb),' max=',max(gorb)
;print, 'max if for max(mdata.visphierr)= ',max(mdata.visphierr)
;;;;;;;;;;;;;
;;;calcul de w_rad et w_tan
fact_plus = .5D*(EXP(-2.D*mdata.visphierr^2)+1D)
fact_moins = .5D*(-EXP(-2.D*mdata.visphierr^2)+1D)
fact_carre = .5D*(-EXP(-mdata.visphierr^2)+1D)^2
cmdata.w_rad = (fact_plus*abs2(mdata.visamperr)+fact_carre*abs2(mdata.visamp))^(-1)
cmdata.w_tan = (fact_moins*(abs2(mdata.visamperr)+abs2(mdata.visamp)))^(-1)

; it is VERY important to protect w_rad and w_tan about infinities if
; we want fmin_op no to be stuck in an indefinite loop!!!!
badvalues=where(~FINITE(cmdata.w_rad),count)
if (count GT 0) then begin
   temp=cmdata.w_rad
   temp[badvalues]=0.0
   cmdata.w_rad=temp
endif
badvalues=where(~FINITE(cmdata.w_tan),count)
if (count GT 0) then begin
   temp=cmdata.w_tan
   temp[badvalues]=0.0
   cmdata.w_tan=temp
endif


return,  cmdata
END
;
;
;filename = '~/idl/BIC/PHASE1/final1.oifits'
;data = wisard_oifits2data(filename)
;;seed = 1
;;ndat = 1000
;; data0 = {vis2:abs(randomn(seed, 15))*10., vis2err:abs(randomn(seed, 15)), clot:randomn(seed,10), $
;;          cloterr:abs(randomn(seed, 10)), freqs_u:abs(randomn(seed, 5))+1., freqs_v:abs(randomn(seed, 5))+1.}
;; data=replicate(data0,ndat)
;
;
;mdata = WISARD_DATA2MDATA(data, OPERATORS = operators)
;cmdata = WISARD_MDATA2CMDATA(mdata, OPERATORS = operators)
;H=wisard_make_h(FREQS_U=reform(cmdata.freqs_u,$
;operators.n_bases*n_elements(mdata)),FREQS_V=reform(cmdata.freqs_v, $
;                                                    operators.n_bases*n_elements(mdata)))
;npix = (size(H))[2]
;sqnpix = sqrt(npix)
;TH = conj(transpose(H))
;dirty_map = reform(TH#reform(cmdata.vis,(size(H))[1]), sqnpix, sqnpix)
;x = randomn(seed, npix)
;achix = reform(H#x, operators.n_bases, n_elements(mdata))
;
;END
