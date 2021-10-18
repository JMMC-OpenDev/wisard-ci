; $Id: abs2.pro,v 1.2 2010-10-25 16:56:49 mvannier Exp $
FUNCTION abs2, x, version=version,  help=help

;+
;NOM :
;   ABS2 - Calcul direct du module carré
;   
;CATEGORIE :
;   Mathematics Routines
;
;SYNTAXE :
;	y=abs2(x [, /VERSION] [, /HELP])
;
;DESCRIPTION :
;	Calcule |x|^2 sans passage inutile par SQRT (interne à abs) puis ^2.
;   Le traitement est cependant un peu plus lourd (tests...) que abs(x)^2.
;   Il n'y a donc gain que pour des tableaux, a partir d'environ 1000 éléments.
;   
;   ARGUMENTS :
;    x       : (entrée) la variable a convertir
;
;   /VERSION : (entrée) affichage de la version avant l'exécution.
;   
;   /HELP    : (entrée) affichage de la syntaxe et sortie du programme.
;
;DIAGNOSTIC D'ERREUR :
;   
;
;VOIR AUSSI :
;	abs, imaginary, conj, real
;
;AUTEUR :
;   $Author: mvannier $
;
;HISTORIQUE :
;   $Log: not supported by cvs2svn $
;   Revision 1.3  2002-10-28 10:43:10+01  cassaing
;   Passage () en [], test vals aléatoire au lieu de 0 => abs2 meilleur que abs^2
;
;   Revision 1.2  2001-01-17 14:57:21+01  cassaing
;   Passage doc en doc_library
;
;   Revision 1.1  1999-04-14 18:23:48+02  cassaing
;   Initial revision
;
;   Revision 1.1  1999-04-14 18:17:49+02  cassaing
;   Initial revision
;
;	V3.0, 14/04/1999 - FC : plus de transtypage entier->float, passage RCS
;	V2.0, 28/11/1996 - F. CASSAING : prise en compte tous les types (sug L.M.)
;	V1.1, 06/11/96      F. Cassaing : test type direct (suggestion L.Mugnier)
;	V1.0, 29/07/96      F. Cassaing : 1er jet
;-

on_error,2
IF keyword_set(version) THEN message, '$Revision: 1.2 $, $Date: 2010-10-25 16:56:49 $', /INFO

IF (n_params() NE 1) OR keyword_set(help) THEN  BEGIN
    message, 'Aide demandée ou syntaxe incorrecte. Documentation :', /INFO
    doc_library, routine_courante()
    retall
ENDIF

; le case suivant est classé par ordre d'utilisation a priori décroissant. 
; mais est-ce que IDL ne reorganise pas ?
CASE size(x, /type) OF
;    6: begin &t=abs(x)& return, t*t &  end
    6: return, float(conj(x)*x) ; complex
;    9: begin &t=abs(x)& return, t*t &  end
    9: return, double(conj(x)*x) ; double complex
    4: return, x^2    ; float
    5: return, x^2    ; double
    1: return, x^2    ; byte
    2: return, x^2    ; integer
    3: return, x^2    ; long integer
    ELSE: message, 'type incorrect' ; struct, chaine, undef...
ENDCASE
END                   ; de la fonction abs2

PRO abs2test, chaine, NB
; pour faire ce test, on rempli le tab avec des nbs au hasard.
; car l'algo de sqrt [a priori complexe] peut etre rapide pour x=0
IF execute('var='+chaine) NE  1 THEN message, 'argument invalide'
var[*] = randomn(seed, n_elements(var))
t0 = systime(1) & FOR i = 1, NB DO x = abs(var)^2
t1 = systime(1) & FOR i = 1, NB DO x = abs2(var)
t2 = systime(1) & FOR i = 1, NB DO x = abs(var) & x*=x
t3 = systime(1) & print, chaine, t1-t0, t2-t1, t3-t2
END                   ; du pro abs2test

; debut du main de test
print, 'Tps execution boucle:    abs^2           abs2       abs&x*=x'
abs2test, 'findgen(1024)', 1000
abs2test, 'findgen(65536L)', 100
abs2test, 'float(3)         ', 10000L
abs2test, 'fltarr(128,128)   ', 10000L
abs2test, 'dblarr(128,128)   ', 10000L
abs2test, 'complexarr(16, 16)', 100
abs2test, 'complexarr(32, 32)', 100
abs2test, 'complexarr(128, 128)', 100
abs2test, 'complexarr(256, 256)', 100
abs2test, 'dcomplexarr(128,128)', 100

END                   ; du main de test
