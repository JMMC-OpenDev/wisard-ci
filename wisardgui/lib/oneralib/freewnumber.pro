; $Id: freewnumber.pro,v 1.2 2010-10-25 16:56:49 mvannier Exp $
; Copyright (c) L. Mugnier, ONERA, 2007.
;
; This software is copyright ONERA.
; The author is Laurent Mugnier.
; E-mail: mugnier at onera.fr
;
; This software is a computer program whose purpose is to find the first
; available window under IDL.
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

FUNCTION FREEWNUMBER, VERSION = version, HELP=help
;+
;NOM :
;   FREEWNUMBER - Retourne le numéro de la première fenêtre disponible.
;
;CATEGORIE :
;   Window Routines.
;
;SYNTAXE :
;   fenetre_dispo = FREEWNUMBER( [/HELP] [, /VERSION])
;
;DESCRIPTION :
;   FREEWNUMBER retourne le numéro de la première fenêtre disponible entre 0
;   et 31 (-1 s'il n'y en a pas). Utile pour créer de nouvelles fenêtres. Ex. :
;   IDL> window, Freewnumber()
;   alloue une nouvelle fenêtre.
;
;DIAGNOSTIC D'ERREUR :
;   FREEWNUMBER retourne -1 s'il n'y a pas de fenêtre disponible.
;
;VOIR AUSSI :
;   DEVICE, window_state (retourne l'état des fenêtres de 0 à 31).
;
;AUTEUR :
;   $Author: mvannier $
;
;HISTORIQUE :
;   $Log: not supported by cvs2svn $
;   Revision 1.1  2007/01/26 11:32:48  mugnier
;   Initial revision
;
;   1.1 = version control and /version keyword.
;   Laurent Mugnier - novembre 1995. Mot-clé help.
;   Laurent Mugnier - septembre 1994.
;
;-

on_error,2
IF keyword_set(version) THEN $
    printf, -2, '% '+routine_courante()+': $Revision: 1.2 $, $Date: 2010-10-25 16:56:49 $'
IF (n_params() NE 0) OR keyword_set(help) THEN $
    message, "Usage : available_window = freewnumber( [/HELP] [, /VERSION])"

windowstate = bytarr(32) &  count = 0
device, window_state=windowstate

windowstate = where((windowstate EQ 0),count);indices des elements nuls
IF (count NE 0) THEN  $
  return, windowstate[0] $                   ;1er indice de fenetre disponible
ELSE $
  return,  -1
END
