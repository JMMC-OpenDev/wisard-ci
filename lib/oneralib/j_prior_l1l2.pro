; $Id: j_prior_l1l2.pro,v 1.2 2010-10-25 16:56:49 mvannier Exp $ 
; Copyright (c) L. Mugnier, ONERA, 2007.
;
; This software is copyright ONERA.
; The author is Laurent Mugnier.
; E-mail: mugnier at onera.fr
;
; This software is a computer program whose purpose is to compute the
; criterion for L1-L2 regularization.
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

FUNCTION J_PRIOR_L1L2, object, SCALE = scale, DELTA = delta,  $
                       GRADIENT_O = gradient_o, $
                       WHITE = white, MEAN_O = mean_o, $
                       VERSION = version, HELP = help
;+
;NOM :
;   J_PRIOR_L1L2 - Crit�re pour un a priori L1-L2.
;   
;CATEGORIE :
;   Mathematics Routines
;
;SYNTAXE :
;   critere = J_PRIOR_L1L2(object [, SCALE=scale] , DELTA=delta $
;    [,/WHITE [,MEAN_O=mean_o]] [, GRADIENT_O=gradient_o] [,/VERSION] [,/HELP])
;
;DESCRIPTION :
;
;   Calcule et renvoie la valeur du crit�re J � minimiser pour un a priori
;   L1-L2 stationnaire sur l'objet recherch�. L1-L2 signifie quadratique pour
;   les petites veleurs (du gradient ou de l'objet), lin�aire pour les valeurs
;   >> delta. Calcul optionnel du gradient de ce crit�re.
;   Formule : J(o) = delta^2 . \sum ( |Do|/delta -ln(1+|Do|/delta) ), o� 
;             |Do| est la norme soit du gradient objet [sqrt(o'_x^2 + o'_y^2)] 
;             soit de l'objet (avec /white) normalis� par scale. 
;   
;   Remarques :
;   1 - le crit�re est convexe, positif, nul pour |Do|=0.
;
;   2 - le crit�re est invariant par changement d'�chelle (o -> k.o) si l'on
;       change �galement l'�chelle de scale (scale -> k.scale). 
;       Il cro�t lin�airement avec le nombre de points de l'objet.
;
;   3 - Lien avec le param�trage plus classique et pas invariant :
;       . Formule implant�e, avec delta et scale :
;         (1) J(o) = delta^2 . \sum( |Do|/[delta.scale] - ln(1+|Do|/[delta.scale]))
;       . Formule plus classique et pas invariante, avec delta' et mu :
;         (2) J(o) = mu.deltaprime^2 . \sum(|Do|/[deltaprime] - ln(1+|Do|/[deltaprime]))
;       o� |Do| est la norme soit du gradient objet [sqrt(o'_x^2 + o'_y^2)] 
;       soit de l'objet (avec /white). 
;       . Lien entre les deux (identification) :
;         delta^2 = mu * deltaprime^2  et  delta*scale=deltaprime
;         => mu = 1./scale^2  &  deltaprime = delta*scale
;       . Pour utiliser la formule (2) il faut donc prendre :
;         => delta = sqrt(mu) * deltaprime  &  scale = 1./sqrt(mu)
;            -----------------------------     -------------------
;
;   ARGUMENTS :
;   object    : (entr�e) objet 2D pour lequel on veut la valeur du crit�re.
;
;   SCALE     : (entr�e) facteur d'�chelle appliqu� � Do ; vaut 1 par d�faut.
;               NB : scale^2 joue le r�le de l'inverse d'un hyper-param�tre.
;
;   DELTA     : (entr�e) seuil sur le gradient objet (ou sur l'objet avec
;               /white) pour le changement de comportement
;               quadratique/lin�aire. Quand delta augmente et tend vers
;               l'infini, le crit�re tend vers un crit�re quadratique.
;               Si scale ~ �cart-type du gradient objet (ou de l'objet avec /white),
;               alors prendre delta ~ 2 � 3.
;
;   MEAN_O    : (entr�e) objet moyen, utile si /WHITE pass�.
;
;   GRADIENT_O: (entr�e/sortie) si d�fini et non nul en entr�e, contient le
;               gradient du crit�re au point object en sortie.
;
;   WHITE    : (entr�e) loi a priori "bruit blanc" c�d pixels objets
;              ind�pendants. Le crit�re L1-L2 est alors non sur les gradients
;              de l'objet mais sur les valeurs de l'objet.
;
;   /VERSION : (entr�e) affichage de la version avant l'ex�cution.
;   
;   /HELP    : (entr�e) affichage de la syntaxe et sortie du programme.
;
;DIAGNOSTIC D'ERREUR :
;
;VOIR AUSSI :
;   J_DATA_LS, J_PRIOR_GAUSS
;
;AUTEUR :
;   $Author: mvannier $
;
;HISTORIQUE :
;   $Log: not supported by cvs2svn $
;   Revision 2.3  2007/10/29 17:03:33  mugnier
;   Ajout� dans la documentation le calcul de delta et scale si l'on veut
;   r�gulariser par la formule plus classique :
;   J(o) = mu.deltaprime^2 . \sum(|Do|/[deltaprime] - ln(1+|Do|/[deltaprime])).
;
;   Revision 2.2  2007/01/26 15:22:42  mugnier
;   Added copyright and CeCILL-C license.
;
;   Revision 2.1  2006/04/10 12:38:42  mugnier
;   Nouveau mot-cl� MEAN_O (objet a priori) utilisable avec /WHITE.
;
;   Revision 2.0  2005/01/06 18:05:48  mugnier
;   Mot-cl� /WHITE pour crit�re L1-L2 sur l'objet et non son gradient
;   (mod�le d'objet bruit blanc c�d � pixels ind�pendants).
;
;   Revision 1.3  2003/06/30 13:08:06  mugnier
;   Routine obsolete STDEV remplac�e par STDDEV.
;
;   Revision 1.2  1998-12-22 14:59:28+01  mugnier
;   Introduction du param�tre d'�chelle SCALE comme hyper-param�tre.
;
;   Revision 1.1  1998-11-12 20:28:02+01  mugnier
;   Initial revision
;
;-

on_error,2
IF keyword_set(version) THEN $
    message, "$Revision: 1.2 $, $Date: 2010-10-25 16:56:49 $", /INFO

IF (n_params() NE 1) OR keyword_set(help) OR NOT keyword_set(delta) THEN  $
    message,  $
    'Usage : critere = J_PRIOR_L1L2(object, [, SCALE=scale], DELTA=delta, ' + $
    '[/WHITE [, MEAN_O=mean_o]]' + $
    '[, GRADIENT_O=gradient_o] [, /VERSION] [, /HELP])'

IF NOT keyword_set(scale) THEN scale = 1.
scale_float = 1.0 * scale

delta_float = 1.0 * delta ; calculs faux si delta entier (par ex. 1000^2=16960)

IF keyword_set(white) THEN BEGIN
    IF NOT keyword_set(mean_o) THEN mean_o = 0.0
    DxO = (object-mean_o)/(delta_float*scale_float) ; avec signe, pour gradient.
;    DyO = 0.
    sqrtDxO2plusDyO2 = abs(DxO)
ENDIF ELSE BEGIN
    DxO = (object - shift(object, 1, 0))/(delta_float*scale_float)
    DyO = (object - shift(object, 0, 1))/(delta_float*scale_float)
    sqrtDxO2plusDyO2 = sqrt(DxO^2+DyO^2)
ENDELSE

critere = (delta_float^2) * total(sqrtDxO2plusDyO2-alog(1.+sqrtDxO2plusDyO2))

IF keyword_set(gradient_o) THEN BEGIN
    IF keyword_set(white) THEN $
        gradient_o = (delta_float/scale_float) * (DxO / (1.+sqrtDxO2plusDyO2)) $
    ELSE BEGIN
        DxO = temporary(DxO) / (1. + sqrtDxO2plusDyO2);economie d'un tableau suppl.
        DyO = temporary(DyO) / (1. + sqrtDxO2plusDyO2);economie d'un tableau suppl.
        gradient_o = (delta_float/scale_float) * $
                     (DxO - shift(DxO, -1, 0) + DyO - shift(DyO, 0, -1))
    ENDELSE
    
ENDIF

return, critere

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

Np = 32
object= rebin(readfits('/home/opc/images/galaxie/m51-64x64.fits'), NP, NP, /SAM)

gradient = 1 & dummy = j_prior_l1l2(object,gradient_o=gradient, scale=10., delta=3.)
gradnum=grad_diff_finies('j_prior_l1l2', object, scale=10., delta=3., epsilon = .1)
; NB : il faut ajuster le epsilon a la main
print, 'Corr�lation entre les 2 gradients :', correlate(gradient,gradnum)
print, 'Rapport d''�chelle :', stddev(gradient)/stddev(gradnum)

;; avec /WHITE :
gradient = 1 & dummy = j_prior_l1l2(object,gradient_o=gradient, scale=10., $
                                    delta=3., /white, mean_o = object/2)
gradnum=grad_diff_finies('j_prior_l1l2', object, scale=10., delta=3., /white, $
             mean_o = object/2, epsilon = .1); NB : il faut ajuster le epsilon a la main
print, 'Corr�lation entre les 2 gradients (avec /white) :', $
       correlate(gradient,gradnum)
print, 'Rapport d''�chelle (avec /white) :', stddev(gradient)/stddev(gradnum)

aff, [gradient, gradnum, 1e2*(gradient - gradnum)]

END


