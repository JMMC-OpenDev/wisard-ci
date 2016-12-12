; $Id: j_prior_gauss.pro,v 1.2 2010-10-25 16:56:49 mvannier Exp $ 
; Copyright (c) L. Mugnier, ONERA, 2007.
;
; This software is copyright ONERA.
; The author is Laurent Mugnier.
; E-mail: mugnier at onera.fr
;
; This software is a computer program whose purpose is to compute the
; criterion for quadratic (gaussian) regularization.
;
; This software is governed by the CeCILL-C license under French law and
; abiding by the rules of distribution of free software. You can use, modify
; and/ or redistribute the software under the terms of the CeCILL-C license as
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
; knowledge of the CeCILL-C license and that you accept its terms.
;
FUNCTION J_PRIOR_GAUSS, object, FT_OBJECT=ft_object, PSD = psd, $
          FT_MEAN_O = ft_mean_o, GRADIENT_O = gradient_o, $
          VERSION = version, HELP = help
;+
;NOM :
;   J_PRIOR_GAUSS - Critère pour un a priori gaussien stationnaire.
;   
;CATEGORIE :
;   Mathematics Routines
;
;SYNTAXE :
;   critere = J_PRIOR_GAUSS([object | FT_OBJECT=ft_object], PSD=psd, $
;      [FT_MEAN_O=ft_mean_o] [, GRADIENT_O=gradient_o] [, /VERSION] [, /HELP])
;
;DESCRIPTION :
;
;   Calcule et renvoie la valeur du critère J à minimiser pour un a priori
;   gaussien stationnaire sur l'objet recherché. C'est l'opposé de la
;   log-probabilité de la loi a priori. Calcul optionnel du gradient de ce
;   critère.
;   
;   ARGUMENTS :
;
;   object    : (entrée) objet 2D pour lequel on veut la valeur du critère. Au
;               lieu de object on peut passer le mot-clé FT_OBJECT.
;
;   FT_OBJECT : (entrée) DFT de l'objet, supposée centrée en (0,0). À passer
;               au lieu de object pour un calcul plus rapide.
;
;   PSD       : (entrée) Densité Spectrale de Puissance de l'objet. Si l'objet
;               est de somme Nphotons, alors PSD(0,0) doit valoir Nphotons^2.
;
;   FT_MEAN_O : (entrée) DFT de l'objet moyen, supposée centrée en (0,0).
;               Vaut 0 par défaut.
;
;   GRADIENT_O: (entrée/sortie) si défini et non nul en entrée, contient le
;               gradient du critère au point object en sortie.
;
;   /VERSION : (entrée) affichage de la version avant l'exécution.
;   
;   /HELP    : (entrée) affichage de la syntaxe et sortie du programme.
;
;DIAGNOSTIC D'ERREUR :
;
;VOIR AUSSI :
;   J_DATA_LS.
;
;AUTEUR :
;   $Author: mvannier $
;
;HISTORIQUE :
;   $Log: not supported by cvs2svn $
;   Revision 1.7  2007/01/26 15:24:25  mugnier
;   Added copyright and CeCILL-C license.
;
;   Revision 1.6  2007/01/24 11:27:51  mugnier
;   On appelle désormais la FFT IDL et plus VFFT car la FFT IDL est devenue bien
;   plus rapide (à partir de IDL version 6.0) et autant que fftw a priori.
;
;   Revision 1.5  2003/06/30 13:06:58  mugnier
;   Routine obsolete STDEV remplacée par STDDEV.
;
;   Revision 1.4  2000-06-09 13:58:25+02  woillez
;   Corrigé bug dim2 pour images de taille >=256.
;
;   Revision 1.3  2000-05-30 15:36:27+02  mugnier
;   Appel à VFFT (Very FFT) en lieu et place de la fonction IDL de base.
;
;   Revision 1.2  1998-12-02 10:44:18+01  mugnier
;   Utilisation de REAL() plutôt que de FLOAT() : préserve le type des entrées.
;
;   Revision 1.1  1998-11-12 20:24:11+01  mugnier
;   Initial revision
;
;-

on_error,2
IF keyword_set(version) THEN $
    message, "$Revision: 1.2 $, $Date: 2010-10-25 16:56:49 $", /INFO

IF (n_params() GT 1) OR keyword_set(help) THEN message,  $
  'Usage : critere = J_PRIOR_GAUSS([object | FT_OBJECT=ft_object], PSD=psd, [FT_MEAN_O=ft_mean_o] [, GRADIENT_O=gradient_o] [, /VERSION] [, /HELP])'

; NORMALISATION DSP : 
; si total(object) = Nph, alors DSP(0,0) doit valoir Nph^2. 
; Explique l'un des dim2 du critère, et le dim2 du gradient.
; L'autre dim2 du critère vient de la TFD d'IDL

IF (n_params() EQ 1) THEN  $
    ft_object = fft(object, -1)  $
ELSE IF NOT keyword_set(ft_object) THEN message, 'ft_object absent'

IF NOT keyword_set(ft_mean_o) THEN ft_mean_o = 0.
IF NOT keyword_set(psd) THEN message, 'PSD absente.'
dim2 = 1. * n_elements(ft_object)

IF keyword_set(gradient_o) THEN  $
    gradient_o = dim2 * real(fft((ft_object - ft_mean_o)/psd, 1))

return, (dim2^2/2.) * total( abs2(ft_object - ft_mean_o)/psd )

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; object supposé défini et 32x32 :
gradient = 1 & dummy = j_prior_gauss(object,gradient_o=gradient, psd = 1.)
gradnum=grad_diff_finies('j_prior_gauss',object, psd = 1.)
print, 'Corrélation entre les 2 gradients :', correlate(gradient,gradnum)
print, 'Rapport d''échelle :', stddev(gradient)/stddev(gradnum)

END
