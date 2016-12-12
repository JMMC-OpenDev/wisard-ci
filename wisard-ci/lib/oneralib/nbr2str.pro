; $Id: nbr2str.pro,v 1.3 2010-11-30 21:26:10 duvert Exp $
; Copyright (c) L. Mugnier, ONERA, 2007.
;
; This software is copyright ONERA.
; The author is Laurent Mugnier.
; E-mail: mugnier at onera.fr
;
; This software is a computer program whose purpose is to convert an array of
; numbers to the corresponding strings.
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
FUNCTION NBR2STR, numbers, FORMAT = format, PRINT = print, $
                  VERSION = version, HELP = help

;+
;NOM :
;   NBR2STR - Conversion de (tableau de) nombre(s) en chaîne de caratères.
;   
;CATEGORIE :
;   String Processing Routines
;
;SYNTAXE :
;   chaine = NBR2STR(numbers [, FORMAT=format] [, /STRING]
;                    [, /VERSION] [, /HELP])
;
;DESCRIPTION :
;   Conversion de (tableau de) nombre(s) en chaîne de caratères.
;   
;   ARGUMENTS :
;   numbers  : nombre ou vecteur de nombres à convertir en chaîne.
;
;   FORMAT   : (entrée) format utilisé par la fonction STRING pour formatter
;              l'entrée. Par défaut, FORMAT='(F25.2)' pour réels et complexes.
;
;   /PRINT   : (entrée) si présent, passé à la fonction STRING pour formatter
;              l'entrée. Sert à forcer un résultat de type string plutôt
;              que vecteur de string quand on passe un vecteur de nombres en
;              entrée. 
;
;   /VERSION : (entrée) affichage de la version avant l'exécution.
;   
;   /HELP    : (entrée) affichage de la syntaxe et sortie du programme.
;
;DIAGNOSTIC D'ERREUR :
;   
;
;VOIR AUSSI :
;   strcompress
;   string
;
;AUTEUR :
;   $Author: duvert $
;
;HISTORIQUE :
;   $Log: not supported by cvs2svn $
;   Revision 1.5  2009/03/18 14:20:28  lmugnier
;   Removed the dependency on CODE.PRO.
;
;   Revision 1.4  2007/01/26 11:28:20  mugnier
;   Added copyright and CeCILL-C license in source.
;
;   Revision 1.3  2006/05/02 15:22:58  mugnier
;   Mot-clé /PRINT, passé à STRING.
;
;   Revision 1.2  2001-05-18 17:27:30+02  mugnier
;   Mot-clé FORMAT pour choisir le nombre de chiffres après la virgule.
;
;   Revision 1.1  1999-11-26 14:31:15+01  mugnier
;   Initial revision
;
;-

on_error,2
IF keyword_set(version) THEN $
    message, '$Revision: 1.3 $, $Date: 2010-11-30 21:26:10 $', /INFO 

IF (n_params() NE 1) OR keyword_set(help) THEN $ 
    message, 'Usage : chaine = NBR2STR(numbers [, /VERSION] [, /HELP])'

sz = size(numbers)
type_entree = sz[sz[0]+1];data type: 4=float, 5=double, 6=complex, 9=double-comp

IF ((type_entree EQ 4) OR (type_entree EQ 5) OR  $;    réel simple ou double
   (type_entree EQ 6) OR (type_entree EQ 9)) $ ;complexe simple ou double
    AND NOT keyword_set(format) THEN format = '(F25.2)' 
; on veut 2 chiffres apres virgule par defaut.
; pour grands nbs utiliser E25.2 (ou qq chose comme g25.7)

IF keyword_set(print) THEN $
    return, strcompress(string(numbers, FORMAT = format, /PRINT)) $
ELSE $
    return, strcompress(string(numbers, FORMAT = format), /remove_all)

END

