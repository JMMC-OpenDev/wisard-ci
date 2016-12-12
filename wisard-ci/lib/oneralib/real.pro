; $Id: real.pro,v 1.2 2010-10-25 16:56:50 mvannier Exp $ 
; Copyright (c) L. Mugnier, ONERA, 2007.
;
; This software is copyright ONERA.
; The author is Laurent Mugnier.
; E-mail: mugnier at onera.fr
;
; This software is a computer program whose purpose is to extract the real
; part of a complex number or array.
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
FUNCTION REAL, z, VERSION = version, HELP = help

;+
;NOM :
;   REAL - partie réelle d'un nombre ou d'un vecteur de nombres complexe(s)
;   
;CATEGORIE :
;   Type Conversion Routines
;
;SYNTAXE :
;   partie_reelle = REAL(z [, /VERSION] [, /HELP])
;
;DESCRIPTION :
;
;   REAL() calcule et renvoie la partie réelle d'un nombre complexe ou d'un
;   vecteur de nombres complexes, *sans* en changer la précision. REAL() a
;   donc le même comportement que IMAGINARY(), contrairement à FLOAT().
;   
;   Une entrée réelle (simple ou double précision) est tolérée. REAL()
;   retourne alors l'entrée non modifiée.
;   
;   ARGUMENTS :
;   z        : (entrée) nombre ou vecteur de nombres complexe(s) dont on veut
;              la partie réelle.
;
;   /VERSION : (entrée) affichage de la version avant l'exécution.
;   
;   /HELP    : (entrée) affichage de la syntaxe et sortie du programme.
;
;DIAGNOSTIC D'ERREUR :
;   
;EXEMPLE :
;
;IDL> help,real(dcomplex(1,2)) -> <Expression>    DOUBLE    =        1.0000000
;IDL> help,real(complex(1,2))  -> <Expression>    FLOAT     =       1.00000

;VOIR AUSSI :
;   
;
;AUTEUR :
;   $Author: mvannier $
;
;HISTORIQUE :
;   $Log: not supported by cvs2svn $
;   Revision 1.2  2007/01/26 11:26:41  mugnier
;   Ajout copyright et licence CeCILL-C.
;
;   Revision 1.1  1998/11/30 15:40:25  mugnier
;   Initial revision
;
;-

on_error,2
IF keyword_set(version) THEN $
    message, "$Revision: 1.2 $, $Date: 2010-10-25 16:56:50 $", /INFO

IF (n_params() NE 1) OR keyword_set(help) THEN $
    message, "Usage : real_part = REAL(z [, /VERSION] [, /HELP])" 

sz = size(z)
type = sz[sz[0]+1]

;debug : print, type

IF (type EQ 6) OR (type EQ 4) THEN $ ; complexe ou réel simple précision
    return, float(z) $
ELSE IF (type EQ 9) OR (type EQ 5) THEN $ ; complexe ou réel double précision
    return, double(z) $
ELSE message, 'Z is neither real nor complex.'

; solution + directe mais + chère : return, IMAGINARY(COMPLEX(0., 1.) * z)

END

