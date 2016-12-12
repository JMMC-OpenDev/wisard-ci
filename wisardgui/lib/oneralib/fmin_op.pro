; $Id: fmin_op.pro,v 1.3 2010-11-30 21:26:10 duvert Exp $
; Copyright (c) L. Mugnier, ONERA, 2007.
;
; This software is copyright ONERA.
; The author is Laurent Mugnier.
; E-mail: mugnier at onera.fr
;
; This software is a computer program whose purpose is to be an IDL frontend
; to the OptimPack package by Eric Thiébaut (CRAL) that is (almost) a
; plug-and-play replacement for previous minimization engines used at ONERA.
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
PRO FMIN_OP, x, fx, FUNC = func, INIT = init, $
             WS = ws, GRADIENT = gradient, $ ; entrees/sorties necessaires
             ERRF = errf, $                  ; sortie utile
             JOB = job, $       ; sortie optionnelle
             ITER = iter, $     ; sortie optionnelle
             REINIT = reinit, $ ; entree optionnelle
             CONV_THRESHOLD = conv_threshold, $ ; entree optionnelle
             ABSOLUTE_CONV = absolute_conv, $   ; entree optionnelle
             ITMAX = itmax, $                   ; entree optionnelle
             _REF_EXTRA = extra, VISU = visu, $ ; entrees optionnelles
             ACTIVE_SET = active_set, $ ; entree pour contrainte de positivite
             XMIN = xmin, XMAX = xmax, $
             MEMORY = memory, $                 ; entree optionnelle
             LIBRARY = library, ALLTHEWAY = alltheway, $
             VERSION = version, HELP = help
;+
;NOM :
;   FMIN_OP - Minimisation de critère par appel à OptimPack (VMLM avec Bornes)
;   
;CATEGORIE :
;   Mathematics Routines
;
;SYNTAXE :
;   PRO FMIN_OP, x [, fx] , FUNC = func, INIT = init, $
;                WS = ws, GRADIENT = gradient, $;entrees/sorties necessaires
;                [ERRF = errf, $]                    ; sortie utile
;                [JOB = job, $]                      ; sortie optionnelle
;                [ITER = iter, $]                    ; sortie optionnelle
;                [REINIT = reinit, $]                ; entree optionnelle
;                [CONV_THRESHOLD = conv_threshold, $]; entree optionnelle
;                [/ABSOLUTE_CONV]                    ; entree optionnelle
;                [ITMAX = itmax, $]                  ; entree optionnelle
;                [_REF_EXTRA = extra][, VISU = visu, $];entrees optionnelles
;                [ACTIVE_SET = active_set, $]        ; entree optionnelle
;                [XMIN = xmin, ][XMAX = xmax, $]     ; entrees optionnelles
;                [MEMORY = memory, $]                ; entree optionnelle
;                [LIBRARY=library, $]                ; entree optionnelle
;                [/VERSION] [, /HELP]
;
;DESCRIPTION :
;   FMIN_OP minimise une fonction de N variables par une méthode de type VMLM
;   avec Bornes. C'est un ``wrapper'' ou frontal des routines OP_INIT,
;   OP_VMLMB_SETUP et OP_VMLMB développées par Eric Thiébaut (package
;   OptimPack) destiné à remplacer FMIN_GPA (syntaxe identique). FMIN_OP
;   doit être appelé dans une boucle (car effectue une itération de
;   minimisation) sauf avec le mot-cle /ALLTHEWAY.
;
;   Terminaison : la minimisation s'arrête dès que l'une des conditions
;   suivantes est remplie :
;   - on a passé /ALLTHEWAY et convergence atteinte à CONV_THRESHOLD près
;     (i.e. JOB GE 3) ;
;   - on n'a pas passé /ALLTHEWAY et fin d'une itération atteinte 
;     (i.e. JOB GE 2) ;
;   - le nombre d'itérations maximum spécifié (ITMAX) est atteint ;
;   - on a passé /ALLTHEWAY et l'utilisateur tape 'Q' (la lettre Q majuscule).
;   
;   ARGUMENTS :
;
;     x             : (entrée/sortie) valeur courante de l'estimée du point
;                     minimum (en 1ère entrée, point de départ de la
;                     minimisation).
;
;    fx             : (entrée/sortie optionnelle) nom de variable, qui
;                     contient en sortie la valeur f(x) du critère au point x.
;                     Utilisé en entrée comme valeur précédente de f(x)
;                     pour calcul de ERRF. Initialisé à
;                     double((machar()).xmax) si absent en entrée.
;
;    FUNC           : (entrée) nom de la fonction qui calcule le critère et
;                     son gradient au point courant x. La fonction FUNC(x)
;                     doit avoir exactement 2 arguments "positionnel", le
;                     premier étant x (en entrée) et le second étant
;                     gradient_de_FUNC_en_x, qui doit être *calculé* par FUNC.
;                     FUNC peut accepter un nombre quelconque de paramètres
;                     auxiliaires d'entrée, passés par mot-clé ; en
;                     effet tous les mots-clés passés à FMIN_OP que cette
;                     routine ne connait pas sont transférés à FUNC via le
;                     mot-clé spécial _REF_EXTRA (cf doc IDL).
;                     L'appel à FUNC doit donc être de cette forme :
;                     Y= FUNC(x, gradient, KEYWORD1=keyword1, KEYWORD2=keyword2,...)
;
;    INIT           : (entrée) à mettre a 1 lors du premier appel puis à 0
;                     pour toute la suite de la minimisation. 
;                     Inutile avec le mot-clé /ALLTHEWAY.
;
;    WS             : (entrée/sortie) nom de variable qui reçoit en sortie et
;                     pour les itérations suivantes le ``WorkSpace'' pour
;                     VMLM-B (cf. op_vmlmb.pro). 
;                     Inutile avec le mot-clé /ALLTHEWAY. 
;
;    GRADIENT       : (entrée/sortie) nom de variable qui reçoit en sortie et
;                     pour les itérations suivantes le gradient courant.
;                     Inutile avec le mot-clé /ALLTHEWAY.
;
;    ERRF           : (sortie optionnelle) utile contient la valeur du test de 
;                     convergence sur FUNC(x), qui est l'évolution du critère 
;                     FUNC(x) entre les 2 dernières itérations, *relativement* 
;                     à la valeur moyenne du critère (sauf avec
;                     /ABSOLUTE_CONV, voir ce mot-clé).
;
;    JOB            : (sortie optionnelle) utile si on veut l'avis de op_vmlmb
;                     sur la convergence. Contient la valeur de JOB rendue
;                     par op_vmlmb (voir ce programme). Doit valoir 2 si
;                     itération bien terminée, 3 si op_vmlmb considère qu'il a
;                     convergé (cf CONV_THRESHOLD), et 4 s'il n'arrive plus à
;                     avancer.
;
;    ITER           : (sortie optionnelle) contient le nombre d'itérations
;                     effectuées. 
;
;    REINIT         : (entrée optionnelle) permet de réinitialiser
;                     l'algorithme depuis le point courant si présent et
;                     non nul (méthode pour 1ère descente = gradient simple). 
;                     Par défaut on ne réinitialise jamais.
;
;    CONV_THRESHOLD : (entrée optionnelle) utile si on veut l'avis de
;                     op_vmlmb sur la convergence. Si présent, utilisé par
;                     op_vmlmb comme seuil sur l'évolution du critère (en 
;                     relatif par défaut, ou en absolu avec /ABSOLUTE_CONV).
;                     Mot-clé indispensable avec /ALLTHEWAY.
;
;    ABSOLUTE_CONV  : (entrée optionnelle) si ce mot-clé est présent et
;                     non nul en entrée, alors ERRF donne l'évolution du
;                     critère FUNC(x) entre les 2 dernières itérations, en
;                     absolu et non en relatif à la valeur du critère.
;
;    ITMAX          : (entrée optionnelle) nombre maximum d'itérations. A
;                     n'utiliser que si l'on ne veut pas minimiser
;                     complètement le critère. 
;
;    VISU           : (entrée optionnelle) si ce mot-clé est présent et non
;                     nul, visu d'infos pendant la minimisation, toutes les
;                     VISU itérations (ainsi que quand on va sortir.
;                     Les infos visualisées sont job, le nb d'itérations,
;                     le nb d'évaluation du critère et de son gradient,
;                     la valeur du critère f(x) et la valeur du test de
;                     convergence (cf mot-clé CONV_THRESHOLD). 
;    
;    ACTIVE_SET     : (entrée/sortie optionnelle) pour demander une
;                     minimisation sous contrainte de positivité ; tableau de
;                     type BYTE de la taille de x qui doit valoir 1B lors du
;                     premier appel, et vaut 1B sur les points actifs ensuite.
;                     Si en plus d'ACTIVE_SET on passe XMIN on obtient
;                     une minimisation sous la contrainte X >= XMIN (au lieu
;                     de X >= 0).
;                     Si en plus d'ACTIVE_SET on passe XMAX on obtient
;                     une minimisation sous la contrainte X <= XMAX.
;                     
;
;    MEMORY         : (entrée optionnelle) détermine la taille mémoire
;                     utilisée pour approximer l'inverse du Hessien par la
;                     méthode "Variable metric with Limited Memory" i.e. BFGS
;                     à mémoire limitée. MEMORY=5 par défaut (et est
;                     proportionnelle au nombre de gradients conservés en
;                     mémoire). 
; 
;   LIBRARY=library : (entrée optionnelle) nom complet de la bibliothèque
;                     OptimPack_IDL${OSTYPE} (avec le chemin mais sans le .so,
;                     comme requis par OP_INIT).
;                     Sur les machines 64 bits, c'est OP_INIT qui ajoutera _64
;                     au nom de librairie à charger, donnant
;                     OptimPack_IDL${OSTYPE}_64.so (cf doc de OP_INIT).
;                     Exemple = OptimPack_IDLlinux sous linux,
;                     OptimPack_IDLsolaris sous solaris.
;                     Ce mot-clé est normalement *inutile* car cette
;                     bibliothèque est trouvée automatiquement dès qu'elle est
;                     dans le !PATH.
;
;   /ALLTHEWAY      : (entrée optionnelle) si ce mot-clé est présent et
;                     non nul, FMIN_OP itère jusqu'à ce que (JOB GE 3)
;                     i.e. fin de convergence au lieu de (JOB GE 2) par défaut
;                     i.e. fin d'une itération.
;                     Très utile si on ne veut rien afficher entre 2 itérations.
;                     Avec /ALLTHEWAY, les mots-clés INIT, WS et GRADIENT ne
;                     sont plus obligatoires.
;                     CONV_THRESHOLD devient indispensable.
;   
;   /VERSION        : (entrée) affichage de la version avant l'exécution.
;                     Si (VERSION GE 2) alors VERSION est également passé à
;                     FUNC pour avoir sa version (ce qui suppose que FUNC
;                     accepte le mot-clé VERSION).
;   
;   /HELP           : (entrée) affichage de la syntaxe et sortie du programme.
;
;DIAGNOSTIC D'ERREUR :
;   
;EXEMPLE :
;
;   active_set = byte(x) * 0B + 1B ; où x est le "guess" initial
;   FMIN_OP, x, criterion_value, FUNC = 'criterion_computation_function', $
;            /ALLTHEWAY, CONV_THRESHOLD = (machar()).eps, $
;            ITMAX = 1000L,  ACTIVE_SET = active_set, $
;            FUNC_KEYWORD1=func_keyword1, ... 
;
;
;VOIR AUSSI :
;   OP_INIT, OP_VMLMB_SETUP, OP_VMLMB et OP_VMLMB_MSG (Eric Thiébaut),
;   WHEREIS (Laurent Mugnier).
;
;AUTEUR :
;   $Author: duvert $
;
;HISTORIQUE :
;   $Log: not supported by cvs2svn $
;   Revision 2.4  2009/03/25 15:16:08  mugnier
;   Possibilite' d'nterruption propre des ite'rations par 'Q'
;   si on a passÃ© /ALLTHEWAY.
;
;   Revision 2.3  2008/10/24 17:50:52  mugnier
;   - FMIN_OP trouve de'sormais la librairie OptimPack quel que soit l'OS, pourvu
;     qu'elle soit dans le !PATH.
;   - itmax n'est plus affiche' s'il n'est pas passe'.
;
;   Revision 2.2  2008/09/11 13:52:25  mugnier
;   Test de l'existence de la variable $OSTYPE (pas dÃ©finie sur une vieille
;   machine sous bash).
;
;   Revision 2.1  2007/10/30 17:28:18  mugnier
;   Correction de la doc et de l'affichage d'infos sur la minimisation.
;
;   Revision 2.0  2007/10/10 12:09:34  mugnier
;   Version 2.0 : - /ALLTHEWAY est stabilisé
;                   (affichage infos minimisation,
;                   correction petits bugs,...)
;                 - mot-clé ABSOLUTE_CONV.
;                 - calcul correct nb itérations.
;
;   Revision 1.17  2007/10/01 14:16:26  mugnier
;   Amélioration de la doc.
;
;   Revision 1.16  2007/04/02 15:12:12  mugnier
;   Amélioration de la doc et exemple.
;
;   Revision 1.15  2007/01/31 15:04:55  mugnier
;   Doc sur infos affichées.
;
;   Revision 1.14  2007/01/29 15:55:37  mugnier
;   VERSION peut désormais être passé à FUNC (voir doc).
;
;   Revision 1.13  2007/01/26 14:29:39  mugnier
;   Mot-clé VERSION passé à WHEREIS.
;
;   Revision 1.12  2007/01/26 11:46:13  mugnier
;   Added copyright and CeCILL-C license in source.
;
;   Revision 1.11  2007/01/23 17:43:11  mugnier
;   Correction d'un bug si le chemin de la bibliothèque OptimPack contient "../"
;
;   Revision 1.10  2007/01/19 17:50:15  mugnier
;   Correction petit bug affichage info (si on n'en veut pas).
;
;   Revision 1.9  2007/01/16 17:34:51  mugnier
;   INIT=1 n'est plus nécessaire lors d'un appel avec le mot-clé /ALLTHEWAY.
;
;   Revision 1.8  2006/10/31 15:12:10  mugnier
;   Nouveau mot-clé ITMAX.
;
;   Revision 1.7  2006/04/27 14:55:14  mugnier
;   Mot-clé /ALLTHEWAY pour minimisation complète au lieu d'une seule itération.
;
;   Revision 1.6  2006/04/27 13:12:08  mugnier
;   Doc sur library mise à jour conjointement à op_init.
;
;   Revision 1.5  2005/06/08 12:56:19  mugnier
;   La bibliothèque recherchée par défaut est désormais OptimPack_IDL${OSTYPE},
;   soit par exemple OptimPack_IDLlinux ou OptimPack_IDLsolaris.
;   Ceci permet de garder au même endroit des bibliothèques pour différents OSs.
;
;   Revision 1.4  2004/02/12 09:36:31  mugnier
;   Nouveau mot-clé REINIT pour éventuellement réinitialiser la direction
;   de descente régulièrement si stagnation constatée. *Devrait* être inutile.
;   Implémentation OK d'après E. Thiébaut mais à tester plus complètement.
;
;   Revision 1.3  2003/10/02 15:17:20  mugnier
;   - FX peut etre désormais utilisé en entrée (cf. doc) ;
;   - mots-clés XMIN et XMAX (bornes autres que 0 et infinité) ;
;   - correction d'un effet de bord sur le mot-clé VISU.
;
;   Revision 1.2  2003/04/22 15:51:46  mugnier
;   Message d'erreur si version d'IDL < 5.4 : OptimPack teste la valeur de
;     !version.memory_bits.
;   Corrigé bug si conv_threshold était un nom de variable sans variable.
;   Le 2ème argument, fx, est désormais optionnel comme avec FMIN_GPA.
;
;   Revision 1.1  2003/04/16 16:49:42  mugnier
;   Initial revision
;
;-

common op_common, libname
on_error,2
IF keyword_set(version) THEN BEGIN 
   version_inside = version
   printf, -2, '% ' + routine_courante() + $
            ': $Revision: 1.3 $, $Date: 2010-11-30 21:26:10 $'
ENDIF ELSE $
    version_inside = 0L

IF (NOT ((n_params() GT 0L) AND (keyword_set(func)) AND $
         ((arg_present(ws) AND arg_present(gradient)) $
          OR keyword_set(alltheway)))) $
    OR keyword_set(help) THEN BEGIN 
    message, 'Help wanted or incorrect syntax. Documentation :', /INFO
    doc_library, 'fmin_op'
    retall
ENDIF

IF (keyword_set(active_set)) THEN BEGIN
    projection = 1B
    IF ((size(active_set, /type) NE 1L) OR $ ; type 1 = byte
        (n_elements(active_set) NE n_elements(x))) THEN message, $
            'si positivité par projection, active_set doit être bytarr ' + $
            'et de taille de x'
    IF (n_elements(xmax) EQ 0L) THEN xmax = (machar(/double)).xmax
    IF (n_elements(xmin) EQ 0L) THEN xmin = 0.0D
ENDIF ELSE projection = 0B

IF keyword_set(visu) THEN $ ; Ne pas modifier `visu' si passé et 0
                visu_long = visu $
ELSE $
                visu_long = 2147483647L ;=2L^31-1: pas de visu par defaut. 

IF keyword_set(alltheway) THEN job_endvalue = 3L ELSE job_endvalue = 2L
IF n_elements(itmax) EQ 0L THEN BEGIN
   itmaxgiven = 0B
   itmax = 2147483647L        ; 2^31-1
ENDIF ELSE itmaxgiven = 1B

COMMON fmin_op, job_op, eval_op, iter_op, fx_op, errf_op, negsetcount_op

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; INITIALISATION BOUCLE

IF keyword_set(init) OR keyword_set(alltheway) THEN BEGIN ;initialisation 
    IF (float(!version.release) LT 5.4) THEN $
        message, 'OP_INIT (donc FMIN_OP) ne marche pas avant IDL 5.4.'
    job_op = 1L
    eval_op = 0L
    iter_op = 0L ; iter_op va de 0 à itmax-1 (au plus)
    IF (n_elements(fx) NE 0L) THEN fx_op = double(fx) $
    ELSE fx_op = double((machar()).xmax)
    
    errf_op = 0.0D
    
    IF NOT keyword_set(library) THEN BEGIN 
       OP_INIT, prefix = prefix, basename=basename, suffix = suffix ; computed by op_init
       defaultlibname = basename+suffix; OS-dependent
       
       chemin = FILE_SEARCH(STRSPLIT(!PATH, PATH_SEP(/SEARCH_PATH), $  
                                     /EXTRACT) + '/' + defaultlibname) 
       IF chemin[0] EQ '' THEN $
           message, 'bibliothèque '+defaultlibname+' introuvable !' $
       ELSE $
           libname = chemin[0] ; name of library with full path (op_common).
    ENDIF ELSE $
        OP_INIT, library ; if library given with full path, but without suffix.
    
    IF keyword_set(memory) THEN memory_long = long(memory) $
        ELSE memory_long = 5L
    IF (n_elements(conv_threshold) EQ 0L) THEN conv_threshold = 0.0D
    ws = OP_VMLMB_SETUP(n_elements(x), memory_long, $
                        frtol=(keyword_set(absolute_conv) ? $
                               0D : double(conv_threshold)), $
                        fatol=(keyword_set(absolute_conv) ? $
                               double(conv_threshold) : 0D))
        
ENDIF

IF keyword_set(reinit) THEN reinit_NOT_done = 1L ELSE  reinit_NOT_done = 0L

last_fx  = fx_op ;  fx_op=valide en entree si job GE 2
continue = 1B

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; BOUCLE PRINCIPALE
REPEAT BEGIN
    IF (job_op EQ 1L) AND keyword_set(reinit) AND reinit_NOT_done THEN BEGIN
        IF keyword_set(memory) THEN memory_long = long(memory) $
            ELSE memory_long = 5L
        IF (n_elements(conv_threshold) EQ 0L) THEN conv_threshold = 0.0D
        ws = OP_VMLMB_SETUP(n_elements(x), memory_long, $
                            frtol=(keyword_set(absolute_conv) ? $
                                   0D : double(conv_threshold)), $
                            fatol=(keyword_set(absolute_conv) ? $
                                   double(conv_threshold) : 0D))
        
        reinit_NOT_done = 0L
    ENDIF

    IF (job_op EQ 1L) THEN BEGIN ;ÉVALUATION DU CRITÈRE ET DU GRADIENT EN X
       
;       last_fx = fx_op       ; =dernier f(x) evalué qd on demande f(x)
                              ; à calculer uniquement en debut/fin d'iter.
                              ; ici ça peut être trop souvent.
        
        ;; PROJECTION DE X à faire avant l'évaluation :
        IF projection THEN BEGIN
           IF (((iter_op+1L) MOD visu_long) EQ 0L) THEN BEGIN
              negative_set = where((x LT xmin) OR (x GT xmax), negsetcount_op)
              print, 'Nb pts < 0 avant projection : ', nbr2str(negsetcount_op)
           ENDIF
            x = (temporary(x) >  xmin) < xmax
        ENDIF
        
        ;; ÉVALUATION DU CRITÈRE et du gradient en x
        IF (version_inside GE 2) THEN BEGIN
           IF keyword_set(extra) THEN $;on veut toujours gradient=>tjrs passé
               fx_op = CALL_FUNCTION(func, x, gradient, _EXTRA=extra, $
                                     VERSION = version_inside) $
           ELSE $
               fx_op = CALL_FUNCTION(func, x, gradient, VERSION=version_inside)
        ENDIF ELSE BEGIN 
           IF keyword_set(extra) THEN $;on veut toujours gradient=>tjrs passé
               fx_op = CALL_FUNCTION(func, x, gradient, _EXTRA=extra) $ ;_STRICT
           ELSE $
               fx_op = CALL_FUNCTION(func, x, gradient)
        ENDELSE

        eval_op += 1L
;        errf_op = ne pas calculer ici
    ; f(x) est valide en sortie ssi job=2 ou 3 (cf op_vmlm.pro).
    ; donc errf est valide en sortie ssi job=2 ou 3 
    ; en effet si job=2 ou 3 en sortie de op_vmlmb on avait job=1 en entrée.
    ENDIF

    ;; APPEL À L'OPTIMISEUR (après projection éventuelle) :
    ;; NB : gradient doit etre present meme si pas actuel (ie meme si job NE 1)
    IF projection THEN BEGIN 
        IF ((job_op GT 1L) OR (eval_op EQ 1L)) THEN BEGIN
            ;; mise à jour de l'active set (job > 1 ou 1ere eval critere)
;            active_set = ((x GT 0.) OR (gradient LT 0.))
            active_set = ((x GT xmin) OR (gradient LT 0.0D)) AND $
                         ((x LT xmax) OR (gradient GT 0.0D))
;; On est sûr qu'on n'a pas ici de valeurs de x <0 parce que si x est modifié
;; par op_vmlmb alors forcément job=1 en sortie donc on applique les bornes
;; juste après.      
            
            IF (((iter_op+1L) MOD visu_long) EQ 0L) THEN BEGIN 
            zero_set = where((x EQ xmin) OR (x EQ xmax), count1)
            zero_and_positivegrad_set = where(active_set EQ 0B, count) 
            print, 'Après proj, nb pts vérifiant Kuhn-Tucker/nuls : ', $
               nbr2str(count)+'/'+nbr2str(count1)+'.'
            ENDIF 
         ENDIF
        
        job_op = OP_VMLMB(x, fx_op, gradient, ws, ACTIVE = active_set)
    ENDIF ELSE $
        job_op = OP_VMLMB(x, fx_op, gradient, ws)
    
    ;; FIN D'UNE ITERATION :
    IF (job_op GE 2L) THEN BEGIN ;; on a fini 1 itération.
       errf_op = (last_fx - fx_op) ; on calcule la convergence 
       IF NOT keyword_set(absolute_conv) THEN $ ; abs. ou relative
           errf_op *= (2.0D/(fx_op+last_fx))
       last_fx = fx_op
       IF ((((iter_op+1L) MOD visu_long) EQ 0L) $ ; toutes les visu_long iterations
           OR $ ; ou, si visu est passé, quand on a convergé / va sortir
           (keyword_set(visu) AND $
            ((job_op GE job_endvalue) OR (iter_op EQ (itmax-1)))) ) THEN $
           printf,-2, $            ; on affiche des infos
                  'FMIN_OP: Job=', nbr2str(job_op), $
                  '; nb iter=', nbr2str(iter_op+1L), $
                  (itmaxgiven ? '/'+nbr2str(itmax) : ''), $
                  '; eval crit=', nbr2str(eval_op), '; f(x)=', $
                  nbr2str(fx_op, format = '(G25.15)'), '; conv=', errf_op
       iter_op += 1L               ; on incremente le compteur d'iter.
    ENDIF
    
    IF keyword_set(alltheway) AND (getenv('DISPLAY') NE '') THEN BEGIN
;get_kbrd doit avoir un terminal (tty) en input
        IF (strupcase(get_kbrd(0)) EQ 'Q') THEN BEGIN ;interruption propre des itérations
            continue = 0B
            message, 'interruption des iterations.', /INFO
            interrupt_keyboard = 1B 
        ENDIF
    ENDIF  

ENDREP UNTIL ((job_op GE job_endvalue) OR (iter_op EQ itmax) OR $
              (NOT continue))
; sortir avec job=2 ou 3, i.e. à la fin d'une iter ou à convergence,
; ou apres itmax iterations (et job GE 2).


;; passage des valeurs en sortie (okazou on les veut)
fx = fx_op                    ; fx rempli en sortie même si vide en entrée 
errf = errf_op
job = job_op
iter = iter_op

END
