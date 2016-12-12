; $Id: angle.pro,v 1.2 2010-10-25 16:56:49 mvannier Exp $ 
; Copyright (c) L. Mugnier, ONERA, 2007.
;
; This software is copyright ONERA.
; The author is Laurent Mugnier.
; E-mail: mugnier at onera.fr
;
; This software is a computer program whose purpose is to compute the phase of
; a complex or complex array.
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
FUNCTION ANGLE, z, VERSION = version, HELP = help

COMPILE_OPT IDL2;idem "DEFINT32, STRICTARR". Cf "Building IDL applis" chap 4. 
;+
;NOM :
;   ANGLE - Angle (i.e., phase) d'un complexe ou vecteur de complexes
;   
;CATEGORIE :
;   Mathematics Routines
;
;SYNTAXE :
;   result = ANGLE(Z [, /VERSION] [, /HELP])
;
;DESCRIPTION :
;
;   ANGLE(Z) donne la phase de Z scalaire ou vecteur complexe (comme ATAN(Z)
;   en IDL 5.4 et antérieures, et comme ATAN(Z,/PHASE) à partir de IDL 5.6).
;   C'est donc l'équivalent d'un ATAN ``backward-compatible'' pour un argument
;   unique complexe. Le nom est dérivé de la fonction Matlab équivalente.
;
;   NB : La sortie de ANGLE est comprise entre -Pi et Pi, comme ATAN(Z,/PHASE).
;   
;   ARGUMENTS :
;
;   Z        : (entrée) scalaire ou vecteur, complexe.
;
;   /VERSION : (entrée) affichage de la version avant l'exécution.
;   
;   /HELP    : (entrée) affichage de la syntaxe et sortie du programme.
;
;DIAGNOSTIC D'ERREUR :
;   
;VOIR AUSSI :
;
;AUTEUR :
;   $Author: mvannier $
;
;HISTORIQUE :
;   $Log: not supported by cvs2svn $
;   Revision 1.3  2007/01/26 11:51:40  mugnier
;   Added copyright and CeCILL-C license in source.
;
;   Revision 1.2  2003/10/21 10:11:43  mugnier
;   Doc ajoutée sur valeurs possibles en sortie : entre -Pi et Pi.
;
;   Revision 1.1  2003-02-21 14:46:43+01  mugnier
;   Initial revision
;
;-

on_error,2
IF keyword_set(version) THEN $
    printf, -2, '% ANGLE: $Revision: 1.2 $, $Date: 2010-10-25 16:56:49 $'

IF (n_params() NE 1) OR keyword_set(help) THEN BEGIN
    message, 'Help wanted or incorrect syntax. Documentation :', /INFO
    doc_library, 'angle'
    retall
ENDIF

sizez = size(z)
typez = sizez[sizez[0]+1]
IF (typez NE 6L) AND (typez NE 9L) THEN $ 
    message, 'The argument Z must be complex.'
                                               
;;ceci ne compile pas en 5.3 & 5.4: % Keyword parameters not allowed in call.
; IF (float(!version.release) GE 5.6) THEN $ 
;     return, atan(z, /phase) $ ; + rapide que atan(imaginary(z), real(z)) 
; ELSE $                        ; IDL v. 5.5 et antérieures :
;     return, atan(imaginary(z), real(z)) 

IF (float(!version.release) LE 5.4) THEN $; /phase est "automatique" avant 5.5
    return, atan(z) $ 
ELSE $
    return, atan(imaginary(z), real(z)) ; + lent que /phase mais
                              ; compile en 5.4 et marche en 5.5.

END

print, 'Angle(complex(1.,1.))  en degrés = 45 ?', $
       angle(complex(1.,1.)) * (180./!PI)
print, 'Angle(dcomplex(1.,1.)) en degrés = 45 ?', $
       angle(dcomplex(1.,1.)) * (180./!DPI)

print, 'Angle(complex(1.,-1.))  en degrés = -45 ?', $
       angle(complex(1.,-1.)) * (180./!PI)
print, 'Angle(dcomplex(1.,-1.)) en degrés = -45 ?', $
       angle(dcomplex(1.,-1.)) * (180./!DPI)

print, 'Angle(complex(-.5,sqrt(3.)/2.))  en degrés = 135 ?', $
       angle(complex(-sqrt(2.)/2.,sqrt(2.)/2.)) * (180./!PI)

END
