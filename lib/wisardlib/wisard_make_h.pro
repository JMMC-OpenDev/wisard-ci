; $Id: wisard_make_h.pro,v 1.12 2010-10-25 16:56:50 mvannier Exp $ 
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

FUNCTION WISARD_MAKE_H, FREQS_U = freqs_u, FREQS_V = freqs_v, $
                        OVERSAMPLING = oversampling, $
                        FOV = fov, NP_MIN = np_min, $
                        ORIGIN = origin, $
                        NP_OUTPUT = np_output, STEP_OUTPUT = step_output, $
                        STEP_INPUT = step_input, VERBOSE = verbose, $
                        VERSION = version, HELP = help

;+
; NAME:
;   WISARD_MAKE_H - Compute matrix H of interferometer's response
;   
; CATEGORY:
;   Array and Image Processing Routines 
;
; CALLING SEQUENCE:
;   H = WISARD_MAKE_H(FREQS_U = freqs_u, FREQS_V = freqs_v,
;                     FOV = fov, NP_MIN = np_min, ORIGIN = origin,
;                     OVERSAMPLING = oversampling,
;                     NP_OUTPUT = np_output, STEP_OUTPUT = step_output,
;                     STEP_INPUT = step_input, 
;                     [, VERBOSE=verbose] [, /VERSION] [, /HELP]
;
; PURPOSE:
;   (to be translated into English)
;   Calcul de la matrice H permettant de passer de l'objet x observé
;   (concaténé dans un grand vecteur) aux visibilités complexes mesurées par
;   un interféromètre : y = Hx.
;   La matrice H est de type et taille DCOMPLEXARR(ND, NP^2) où
;   ND = n_elements(freqs_u).
;   L'élément (l, (j+NP.k)) est calculé ainsi :
;
;   H(l, (j+NP.k)) = d^2 . exp[-2i.Pi.d. (u(l).j + v(l).k) ],
;
;   où d est le pas d'échantillonnage spatial de l'objet.
;   
; POSITIONAL PARAMETERS:
;   None
;
; KEYWORD PARAMETERS:
;
;   FREQS_U  : (entrée) vecteur des abcisses des fréquences spatiales 2D
;              mesurées. 
;
;   FREQS_V  : (entrée) vecteur des ordonnées des fréquences spatiales 2D
;              mesurées. 
;
;   FOV      : (entrée) champ souhaité pour l'objet. 
;              Si absent on prend un champ égal à 1/min(abs(freqs)).
;
;   NP_MIN   : (entrée) Nombre de Points minimum souhaité pour l'objet.
;              Par convention si (NP_MIN LT 0) on prend exactement
;              abs(NP_MIN) points dans le champ FOV.
;
;   ORIGIN   : (entrée) origine plan objet, milieu par défaut (NP_OUTPUT/2).
;
;   OVERSAMPLING: (entrée) autre manière de forcer un suréchantillonnage du
;                 champ FOV : on augmente NP_Nyquist d'un facteur
;                 OVERSAMPLING, où NP_Nyquist est le nb de points requis pour
;                 bien échantillonner le FOV vue la fréquence max des données.
;
;   NP_OUTPUT   : (sortie) Nombre de Points dans le champ. Calculé comme
;                 max(NP_Nyquist, NP_MIN).
;                 La matrice H a donc NP_OUTPUT^2 colonnes.
;
;   STEP_OUTPUT : (sortie) pas d'échantillonnage "d" de la formule ci-dessus.
;                 Calculé comme FOV/NP_OUTPUT
;
;   STEP_INPUT  : (entrée) debug seulement. Pas recommandé.
;
;   VERBOSE     : (entrée) mode bavard (par défaut).
;
;   /VERSION    : (entrée) affichage de la version avant l'exécution.
;   
;   /HELP       : (entrée) affichage de la syntaxe et sortie du programme.
;
; AUTHORS:
; Serge Meimon and Laurent Mugnier (ONERA). 
; Contact: lastname at onera.fr
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
;   freqs = [[1D, 4D], [0D, 3D]] ; ratio 5 between abs([1,0] et [4,3]), so NP_min=10
;   freqs_u = freqs[*, 0]
;   freqs_v = freqs[*, 1]

;   H = WISARD_MAKE_H(FREQS_U = freqs_u, FREQS_V = freqs_v, $
;                     NP_OUTPUT = np_output, STEP_output = step)
;	For more details, see the "main" of test at the end of the routine code. 
;
; SEE ALSO:
;   All the WISARD_*.pro files, which are part of WISARD and covered by the
;   same copyright and license.
;   
; HISTORY:
;   $Log: not supported by cvs2svn $
;   Revision 1.10  2008/01/08 16:41:23  mugnier
;   Fixed typo.
;
;   Revision 1.9  2007/10/31 12:46:35  mugnier
;   Added 3rd reference.
;
;   Revision 1.8  2007/01/30 17:13:10  mugnier
;   Corrected mini-bug in the display of the FOV.
;
;   Revision 1.7  2007/01/29 13:13:38  mugnier
;   Coorected small bugs in display of STEP.
;
;   Revision 1.6  2007/01/26 14:54:19  mugnier
;   Added copyright and CeCILL-B license in source.
;
;   Revision 1.5  2006/10/31 13:52:44  mugnier
;   Nouveau mot-clé ORIGIN.
;
;   Revision 1.4  2006/10/31 10:09:51  mugnier
;   Documentation + mode bavard par défaut.
;
;-

on_error,2
IF keyword_set(version) THEN $
    printf, -2, '% '+routine_courante(VERSION = version) + $
            ': $Revision: 1.12 $, $Date: 2010-10-25 16:56:50 $'

IF (n_params() NE 0) OR keyword_set(help) OR $
    NOT (keyword_set(freqs_u) AND keyword_set(freqs_v))THEN BEGIN 
    message, 'Help required or incorrect syntax. Documentation:', /INFO
    doc_library, routine_courante()
    retall
ENDIF

IF (n_elements(verbose) EQ 0L) THEN verbose = 1B

abs_freqs = sqrt((freqs_u)^2 + (freqs_v)^2)
max_abs_freqs = max(abs_freqs)
min_abs_freqs = min(abs_freqs)


;; FOV est soit passé en arg, soit 1/min(abs(freqs)) :
unemas = 1d-3*(!DPi/180D)/3600D ; une mas en rd

IF keyword_set(fov) THEN $
    fov_inside = fov $
ELSE BEGIN
   IF (float(min_abs_freqs) EQ 0.0) THEN $
       message, 'freqs must not contain the 0 frequency'
   fov_inside = 1D / min_abs_freqs
   IF (verbose) THEN BEGIN
      print, 'FOV = 1/min(abs(freqs)) = ', $
              nbr2str(fov_inside, format = '(G25.7)'), ' = ', $
              nbr2str(fov_inside/unemas), ' mas (if freqs is in rd^{-1}).'
   ENDIF 
   ENDELSE

;; STEP et NP : valeurs MINI données par Nyquist :
IF NOT keyword_set(oversampling) THEN oversampling = 1D
step_nyquist = 1D / (2D * max_abs_freqs * oversampling) ; en radians
np_nyquist   = ceil(fov_inside / step_nyquist)
IF (verbose) THEN $
    print, 'NP min for Nyquist criterion = ', nbr2str(np_nyquist), $
            '. STEP max for Nyquist criterion = ', nbr2str(step_nyquist, $
                                                           format = '(G25.7)')


;; NP final = max(NP nyquist, NP entré) :
IF keyword_set(np_min) THEN BEGIN
   IF np_min GE 0 THEN $      
      np_output = max([np_nyquist, ceil(np_min)]) $ ; long
   ELSE np_output = ceil(abs(np_min)) ; long
; convention : si NP_min<0 alors c'est le NP voulu exactement 
ENDIF ELSE $
    np_output = np_nyquist ; long
IF (n_elements(origin) EQ 0L) THEN origin = np_output/2L

;; STEP final = FOV / NP :
step_output = fov_inside / np_output


IF keyword_set(step_input) THEN $
; uniquement pour debug. Normalement step fixé par FOV etc.
    step_output = double(step_input)

IF (verbose) THEN $
    print, 'NP is finally ', nbr2str(np_output), $
            '. STEP is finally ', nbr2str(step_output, format = '(G25.7)')


nd = n_elements(freqs_u)      ; N. Data

indicerapide = double((lindgen(np_output^2) MOD np_output) - origin)
indicelent = double((lindgen(np_output^2) / np_output) - origin)
; following takes too much memory:
;H = step_output^2 * exp((-2D*!DI*!DPI*step_output) * $
;                        (freqs_u # indicerapide + freqs_v # indicelent))
; this one is a bit longer but saves memory
H=0
H = freqs_u # indicerapide
H = TEMPORARY(H)+freqs_v # indicelent
H = TEMPORARY(H)*(-2D*!DI*!DPI*step_output)
H = exp(TEMPORARY(H))
;H = TEMPORARY(H)*step_output^2 ; removed, and inverse operation in
;wisard.pro, since not used 
; H est de type dcomplexarr, de taille (nd, np_output^2)
return, H

END

;; UNITARY TEST:
; .com wisard_make_h
;NP = 8
;x = dblarr(np, np) & x[0, 0] = 1D
;
;freqs = [[1D, 4D], [0D, 3D]] ; ratio 5 between abs([1,0] et [4,3]), so NP_min=10
;freqs_u = freqs[*, 0]
;freqs_v = freqs[*, 1]
;
;H = WISARD_MAKE_H(FREQS_U = freqs_u, FREQS_V = freqs_v, $
;                  NP_OUTPUT = np_output, STEP_output = step)
;H0 = WISARD_MAKE_H(FREQS_U = freqs_u, FREQS_V = freqs_v, origin = 0)
;Hdefault = WISARD_MAKE_H(FREQS_U = freqs_u, FREQS_V = freqs_v)
;Hnsur2 = WISARD_MAKE_H(FREQS_U = freqs_u, FREQS_V = freqs_v, origin = 5)
;
;H = WISARD_MAKE_H(FREQS_U = freqs_u, FREQS_V = freqs_v, OVER = 1D, NP_MIN = 16, $
;                  NP_OUTPUT = np_output, STEP_output = step, /verbose)
;
;;; Grille frequentielle :
;freqs_u = (lindgen(np^2) MOD np) + 1L
;freqs_v = lindgen(np^2) / np - np/2
;
;H = WISARD_MAKE_H(FREQS_U = freqs_u, FREQS_V = freqs_v, OVER = 1D, NP_MIN = 8, NP_OUTPUT = $
;                  np_output, STEP_output = step)
;
;info, transpose(conj(h))#h ; COMPLEX 324x324
;aff, transpose(conj(h))#h ;, /reim
;
;info, h#transpose(conj(h)) ; DCOMPLEX 64x64, quasi réel, quasi-diag
;aff, h#transpose(conj(h)), /reim
;
;;; Grille frequentielle :
;np = 16L
;freqs_u = (lindgen(np^2) MOD np) - np/2
;freqs_v = lindgen(np^2) / np - np/2
;H = WISARD_MAKE_H(FREQS_U = freqs_u, FREQS_V = freqs_v, OVER = 1D, NP_MIN = $
;                  -np, FOV = 1, NP_OUTPUT = $
;                  np_output, STEP_input = 1D/double(np))
;print, 'NP=', np_output
;print, 'pas=', step, '. 1/pas=', 1D/step
;
;HtH = transpose(conj(h))#h
;info, HtH    ; DCOMPLEX 64x64, quasi réel diag, total 1.0
;aff, HtH
;
;HHt = h#transpose(conj(h))
;info, HHt    ; DCOMPLEX 64x64, quasi réel diag, total 1.0
;aff, HHt
;
;END
