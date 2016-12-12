; $Id: wisard_set_regul.pro,v 2.3 2010-10-25 16:56:50 mvannier Exp $ 
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
;   edited by Jérôme Idier, ISTE, London, 2008.
;


;+
; NAME:
;	WISARD_SET_REGUL - Structure containing regularization information
;
; CATEGORY:
;	Array and Image Processing Routines 
;
; CALLING SEQUENCE:
;	PRIOR = WISARD_SET_REGUL(PSD = psd [, MEAN_O = mean_o], $
;                            DELTA = delta, SCALE = scale [, WHITE = white] $
;                            FWHM = fwhm, MU_SUPPORT = mu_support, $
;                            NP = NP $
;                            [, POSITIVITY = positivity] $
;                            [, SQUARENP = squareNP]
;
; PURPOSE:
;   This function computes and returns an anonymous structure PRIOR containing
;   information on the regularization for use in an image reconstruction
;   routine. 
;   On exit, PRIOR contains the following fields:
;   - type: type of regularization (string). Currently 'l1l2' or 'psd' or
;           'softsupport' or 'none'.
;   - positivity: 0B/1B depending on whether the positivity constraint is to be
;     used or not.
;   - NP: width (Number of Points) of the square object to be reconstructed.
;   - squareNP: square of NP.
;   - some additional fields that depend on `type' and are described below:
;
;   If keyword DELTA (and keyword SCALE) are set on input then type='l1l2'
;   and the additional fields of PRIOR are SCALE, DELTA, WHITE, MEAN_O.
;
;   If keyword PSD is set on input then type='psd' and the additional fields
;   of PRIOR are PSD, FT_MEAN_O.
;
;   If keyword MU_SUPPORT is set on input then type='softsupport' and the
;   additional fields of PRIOR are MU_SUPPORT and MEAN_O.
;
;   These three regularizations are mutually exclusive.
;
; POSITIONAL PARAMETERS:
;	None.
;	
; KEYWORD PARAMETERS:
;
;   PSD        : (input) PSD for PSD (quadratic) regularization, centered on
;                NP/2,NP/2. On ouput, PRIOR.PSD contains PSD centered on (0,0)
;                for use by, e.g., J_PRIOR_GAUSS
;
;   MEAN_O     : (optional input) mean object for quadratic or white L1-L2
;                regularization or soft support for use by, e.g.,
;                J_PRIOR_GAUSS or J_PRIOR_L1L2 or J_PRIOR_SUPPORT.  
;                - For the 2 first regularizations:
;                  If absent, the mean object is taken as a flat object of unit
;                  sum. If you want a null mean object, set MEAN_O explicitly to
;                  0 or to a zero-valued 2D array.
;                  MEAN_O is cast to double precision, normalized to a
;                  unit sum (if not zero) and resampled on NPxNP points. 
;                  On ouput, for PSD regularization, PRIOR.FT_MEAN_O contains
;                  the Fourier transform of MEAN_O centered on (0,0).
;                  For L1-L2 regularization PRIOR.MEAN_O contains MEAN_O.
;                - For the soft support regularization: MEAN_O is used iff
;                  FWHM is not given. it must be a 2D map with stricly
;                  positive values.
;
;   DELTA      : (input) threshold parameter for L1-L2 regularization, for use
;                by, e.g., J_PRIOR_L1L2.
;
;   SCALE      : (input) scale parameter for L1-L2 regularization, for use
;                by, e.g., J_PRIOR_L1L2.
;
;   WHITE      : (optional input) flag telling whether the L1-L2
;                regularization is spike-reserving (independent pixels i.e.,
;                WHITE=1) or edge-preserving (correlated pixels i.e., WHITE=0).
;                WHITE=0 by default. 
;
;   MU_SUPPORT : (input) regularization weight for soft support
;                regularization, for use by J_PRIOR_SUPPORT.
;
;   FWHM       : (optional input) FWHM in pixels for soft support
;                regularization, for use by J_PRIOR_SUPPORT.
;                If given, a Lorentzian is computed as prior object and MEAN_O
;                is discarded even if present. If not given, MEAN_O must be
;                given and is used as the prior object.
;
;   NP         : (input) width of square object to be reconstructed.
;
;   POSITIVITY : (input) flag telling whether POSITIVITY is to be used
;                (1) as a constraint in the reconstruction or not (0).
;                POSITIVITY is set to 1 by default i.e. if absent.
;
;   /VERSION: (input) prints version number before execution.
;   
;   /HELP:    (input) prints the documentation and exits.
;
; AUTHORS:
;   Serge Meimon and Laurent Mugnier (ONERA). 
;   Contact: lastname at onera.fr
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
;   edited by Jérôme Idier, ISTE, London, 2008.
;
; EXAMPLE:
;	gaussianprior = WISARD_SET_REGUL(PSD = psd)
;	spikepreservingprior = WISARD_SET_REGUL(DELTA=0.1, scale=1.0, white=1B)
;	supportprior = WISARD_SET_REGUL(FWHM=10.0, MU_SUPPORT=10.)
;
; SEE ALSO:
;   J_PRIOR_GAUSS, J_PRIOR_L1L2, J_PRIOR_SUPPORT, FMIN_OP.
;   Also see all the WISARD_*.pro files, which are part of WISARD and covered
;   by the same copyright and license.
;
; HISTORY:
;   $Log: not supported by cvs2svn $
;   Revision 2.1  2008/09/19 11:52:03  mugnier
;   Replaced call to 'distc()' by 'dist()', available in any IDL distribution.
;
;   Revision 2.0  2008/03/28 15:29:25  mugnier
;   The "soft support" regularization is now supported
;   (compatible with J_PRIOR_SUPPORT version>=1.2).
;
;   Revision 1.9  2008/01/08 16:40:24  mugnier
;   Fixed typo.
;
;   Revision 1.8  2007/10/31 12:43:48  mugnier
;   Added 3rd reference.
;
;   Revision 1.7  2007/01/29 15:18:25  mugnier
;   MEAN_O = cst (of unit sum) by default but MEAN_O=0 is now possible.
;
;   Revision 1.6  2007/01/25 18:43:49  mugnier
;   Corrected bug if MEAN_O is set and already of the right size.
;
;   Revision 1.5  2007/01/24 17:23:36  mugnier
;   The input PSD is re-sampled on an NPxNP grid and thresholded to strictly
;   positive values if is is not so.
;
;   Revision 1.4  2007/01/18 15:14:56  mugnier
;   Cleanup of code (removed keywords SUPPORT and ACTIVE_SET, cast of MEAN_O in
;   double and POSITIVITY in byte, etc.)
;   +
;   Documentaion.
;
;   Revision 1.3  2007/01/16 17:43:38  mugnier
;   Ajout du champ 'positivity' dans prior.
;
;   Revision 1.2  2007/01/09 11:20:54  meimon
;   ??
;
;   Revision 1.1  2006/11/09 17:57:17  meimon
;   Initial revision
;
;-




FUNCTION WISARD_SET_REGUL, PSD = psd, MEAN_O = mean_o, $
                           DELTA = delta, SCALE = scale, WHITE = white, $
                           MU_SUPPORT = mu_support, FWHM = fwhm, $
                           NP = NP, $
                           POSITIVITY = positivity, $
                           TOTVAR = totvar, $
                           VERSION = version, HELP = help

on_error,2
IF keyword_set(version) THEN $
    printf, -2, '% '+routine_courante()+': $Revision: 2.3 $, $Date: 2010-10-25 16:56:50 $'

IF (n_params() NE 0) OR keyword_set(help) THEN BEGIN 
    message, 'Help required or incorrect syntax. Documentation:', /INFO
    doc_library, routine_courante()
    retall
ENDIF


; POSITIVITÉ PRÉSENTE PAR DÉFAUT :
IF (n_elements(positivity) EQ 0L) THEN $
    positivity_inside = 1B $
ELSE $
    positivity_inside = byte(positivity)

squareNP = long(NP)^2

; OBJET MOYEN :
IF keyword_set(mu_support) THEN BEGIN ; soft support regularization
   IF keyword_set(fwhm) THEN BEGIN ; if FWHM present, discard potential MEAN_O
      distance = double(shift(dist(NP), NP/2L, NP/2L))
      mean_o_inside = 1D / (1D + (2D*distance/fwhm)^2) ;Lorentzian shape
   ENDIF ELSE BEGIN
      IF (n_elements(mean_o) EQ 0L) THEN $
          message, 'For soft support regularization, ' + $
                   'MEAN_O must be given if FWHM is not'
      sizeofmeano = size(mean_o)
      IF ((sizeofmeano[0] NE 2) OR (min(mean_o) LE 0)) THEN BEGIN
         message, 'For soft support regularization, MEAN_O must be 2-D and > 0'
      ENDIF ELSE IF (n_elements(mean_o) NE squareNP) THEN BEGIN ;array, wrong size:
         mean_o_inside = congrid(double(mean_o), NP, NP) > 0 ; resampling may give
                              ; negative values, so threshold result to 0.
      ENDIF ELSE mean_o_inside = double(mean_o) ;if correct size 2D map for mean_o 
   ENDELSE 
ENDIF ELSE BEGIN              ; PSD or L1-L2 regularization
    IF (n_elements(mean_o) EQ 0L) THEN BEGIN ; if absent keyword
       mean_o_inside = dblarr(NP, NP) + 1D/double(squareNP) 
    ENDIF ELSE IF (n_elements(mean_o) EQ 1L) THEN BEGIN ; if scalar keyword
       mean_o_inside = dblarr(NP, NP) + mean_o
    ENDIF ELSE IF (n_elements(mean_o) NE squareNP) THEN BEGIN;if array, wrong size:
       mean_o_inside = congrid(double(mean_o), NP, NP) > 0 ; resampling may give
                                  ; some negative values, so threshold result to 0.
    ENDIF ELSE mean_o_inside = double(mean_o);if correct size 2D map for mean_o 
ENDELSE

total_mean_o = total(mean_o_inside)
IF (total_mean_o NE 0D) THEN $
    mean_o_inside = mean_o_inside / total_mean_o ; unit flux

ft_mean_o = fft(mean_o_inside, -1); pour regul DSP
   
; CRITERE DE REGULARISATION DSP OU L1-L2 (éventuellement blanc)
; on utilise PSD comme hyper en L2 et SCALE en L1-L2.
; la présence ou non de DELTA détermine la régularisation L1-L2 ou L2
; info mise dans la structure "prior"

IF (keyword_set(psd) AND keyword_set(delta)) THEN message,  $
    'PSD (L2) or [SCALE+]DELTA (L1L2), you must choose (at most) one.'
IF (keyword_set(psd) AND keyword_set(mu_support)) THEN message,  $
    'PSD (L2) or MU_SUPPORT (soft support), you must choose (at most) one.'
IF (keyword_set(delta) AND keyword_set(mu_support)) THEN message,  $
    '[SCALE+]DELTA (L1-L2) or MU_SUPPORT (soft support), you must choose (at most) one.'

IF keyword_set(delta) THEN BEGIN ; regularisation L1-L2
   IF (delta GT 100.) THEN message, /INFO, $
        'WARNING: DELTA > 100. Useless (100 is enough) and may give minimization problems!'
   IF (n_elements(scale) EQ 0L) THEN $
       message, 'for an L1-L2 regularization you must input delta and scale'
   IF (n_elements(white) EQ 0L) THEN white = 0L
   
   prior = {type: 'l1l2', $
            scale: scale, delta: delta, white: white, mean_o: mean_o_inside, $
            positivity: (positivity_inside NE 0), NP: NP, squareNP: squareNP }
   
ENDIF ELSE IF keyword_set(psd) THEN BEGIN ; regularization by PSD
   psd_inside = congrid(double(psd), NP, NP) > double((machar(/double)).eps)
   max_psd = max(psd_inside, argmax_psd)
   dummy = where(psd_inside EQ max_psd, count)
   IF ((count EQ 1L) AND (argmax_psd NE NP*NP/2+NP/2)) THEN  $
       message, /INFO, 'The PSD does not seem to be centered.'
   dynamique_psd = max_psd/min(psd_inside)
   IF (dynamique_psd GT 1e6) THEN message, /INFO, $
       'WARNING: very large dynamic range for PSD. Minimization problem may occur!'
    
   prior = {type: 'psd', $
            psd: eclat(psd_inside), $ ; center PSD on (0,0) for j_prior_gauss
            ft_mean_o:ft_mean_o, $
            positivity: (positivity_inside NE 0), NP: NP, squareNP: squareNP } 
ENDIF ELSE IF keyword_set(mu_support) THEN BEGIN ; regularization by soft support
   
   prior = {type: 'softsupport', mu_support: mu_support, $
            mean_o: mean_o_inside, $  ;fwhm: fwhm inutile apres calcul mean_o 
            positivity: (positivity_inside NE 0), NP: NP, squareNP: squareNP } 
ENDIF ELSE IF keyword_set(totvar) THEN BEGIN ; regularization by totvar
   
   prior = {type: 'totvar', positivity: (positivity_inside NE 0), NP: NP, squareNP: squareNP } 
ENDIF ELSE $
    prior = {type: 'none', $
             positivity: (positivity_inside NE 0), NP: NP, squareNP: squareNP }


return, prior

END
    
