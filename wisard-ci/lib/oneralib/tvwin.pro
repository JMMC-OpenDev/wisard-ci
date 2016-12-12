; $Id: tvwin.pro,v 1.2 2010-10-25 16:56:50 mvannier Exp $
pro tvwin,array,win,TITLE=title,XPOS=xpos,YPOS=ypos,FREE=free,SAMPLE=sample $
          ,NOSAMPLE=nosample,ZOOM=zoom,SIZE=size, $
          NOWIN=NOWIN,POSINWIN=posinwin ,REIM=reim,HELP=help,VERSION=version, $
          PROFILE = profile, RDPIX = rdpix  

;+
;PROCEDURE
;	TVWIN - affichage d'une image dans une fenetre de taille adaptee, zoom en option
;
;CATEGORIE : Image Display Routines
;
;SYNTAXE :
;  TVWIN,array[,win,TITLE=title,XPOS=xpos,YPOS=ypos,FREE=free,SAMPLE=sample,NOSAMPLE=nosample,ZOOM=zoom,SIZE=size,NOWIN=NOWIN,POSINWIN=posinwin,REIM=reim,PROFILE = profile, RDPIX = rdpix][, /VERSION] [, /HELP] 
;
;DESCRIPTION
;	TVWIN routine d affichage dans fenetre adaptee
;
;ARGUMENTS
;  ARRAY - Tableau a deux dimensions
;  WIN   - numero de fenetre (par defaut fenetre courante)
;  TITLE, XPOS, YPOS, /FREE: memes arguments que pour commande WINDOW
;  ZOOM  - facteur de zoom en x ou/et y;
;          si ZOOM=-1 -> zoom automatique: petit cote au plus pres de 256 ou
;                                          SIZE si specifie 
;  SIZE  - taille pour le zoom automatique (par defaut 256) 
;  /SAMPLE   - macro-pixels pour le zoom si present (c'est l'option par
;              defaut quand on magnifie)
;  /NOSAMPLE - interpolation lineaire pour le zoom si present (par defaut si
;              on de-magnifie)
;  /NOWIN    - use current window instead of opening new one
;  POSINWIN - position of tvscl in the window (rather interesting if NOWIN set)
;  REIM     - affiche Re(array)/Im(array) si array est un tableau complexe
;              a la meme echelle   si REIM=1
;              en pleine dynamique si REIM=2
;              (par defaut on affiche MODULE/PHASE)
;  PROFILES  - lance profiles dans la foulee
;  RDPIX     - lance rdpix dans la foulee
;  /VERSION : (entrée) affichage de la version avant l'exécution.
;  /HELP    : (entrée) affichage de la syntaxe et sortie du programme.
;
;SEE ALSO
;         WINDOW (idl), FREEWNUMBER (LaMu), TVSCL (idl), BYTSCL (idl), AFF (jmc)
;
;AUTEUR :
;   $Author: mvannier $
;
;HISTORIQUE :
;   $Log: not supported by cvs2svn $
;   Revision 1.7  2001-06-11 17:08:35+02  conan
;   mise a jour des titres par defauts et de la gestion du zoom automatique
;   (utilisation de congrid pour zoom non entier) sur des idees de Laurent Mumu
;
;   Revision 1.6  2001-03-16 16:48:24+01  conan
;   Modif Laurent pour avoir les numeros de fenetres meme en precisant TITLE
;
;   Revision 1.5  1998-05-15 10:59:04+02  conan
;   2 modifs:
;     si on donne SIZE alors ZOOM=-1 par defaut
;     si /FREE on memorise la fenetre courante on y repointe apres affichage
;
;   Revision 1.4  1998-05-13 15:10:50+02  conan
;   ajout de l'option SIZE pour un zoom automatique a taille variable
;
;   Revision 1.3  1998-03-03 16:50:24+01  conan
;   Ajout de la gestion des complexes double precision
;
;   Revision 1.2  1998-03-03 15:56:46+01  conan
;   Mise a jour de la documentation (champs geres par CVS, help, version...)
;
;  17 Fev.  1997 - a la demande de CD&FCh, je change l'option zoom=-1:
;                  zoom identique en x et y, la plus grande dimension est
;                  ramenee au plus proche de 256.
;  30 April 1996 - (1er anniversaire de la fonction !!!) Ajout de la gestion
;                  des tableaux complexes, et en particulier option REIM.
;  08 Janua 1996 - ajout des options NOSAMPLE (les macro-pixels devenant
;						 l'option par defaut), NOWIN, POSINWIN.
;  20 April 1995 - JM CONAN        
;-

on_error,2
IF keyword_set(version) THEN $
    message, "$Revision: 1.2 $, $Date: 2010-10-25 16:56:50 $", /INFO         ;(rempli par CVS ou RCS)

if (n_params() EQ 0) OR keyword_set(help) then message, $
"Usage -TVWIN,array,win,[TITLE,XPOS,YPOS,FREE,SAMPLE,NOSAMPLE,ZOOM,SIZE,NOWIN,POSINWIN,REIM,/VERSION, /HELP]"
 
sz=size(array)
if ( (sz[0] ne 2) ) then $
message,"array must be 2-D"

;-gestion des tableaux complexes-
;type 6 (complexes simple) type 9 (complexes double precision): 
IF (sz[3] EQ 6) OR (sz[3] EQ 9) THEN datacomplex = 1B ELSE datacomplex = 0B

if ( keyword_set(reim) and NOT datacomplex ) then $
message,"Keyword REIM valid only with complex arrays"

if keyword_set(reim) then if ( reim ne 1 ) and ( reim ne 2) then $
message,"Keyword REIM must be equal to 1 or 2 ; other values are not valid"

if (keyword_set(win) and keyword_set(free)) then $
message,"Keywords WIN and FREE can't be used at the same time"

if (keyword_set(nowin) and keyword_set(free)) then $
message,"Keywords NOWIN and FREE can't be used at the same time"

if (keyword_set(nowin) and (keyword_set(title) or (n_elements(xpos) NE 0) or (n_elements(ypos) NE 0))) then $
message,"Keywords TITLE , XPOS and YPOS have no effect if NOWIN is set",/info

if (keyword_set(nowin) and (n_params() eq 2)) then $
wset,win
;message,"Keyword NOWIN is not compatible with a window number"

IF (keyword_set(size) AND NOT keyword_set(zoom)) THEN $
    zoom= -1

IF (keyword_set(sample) AND keyword_set(nosample)) THEN $
message,"Keywords SAMPLE and NOSAMPLE can't be used at the same time"

IF ((keyword_set(sample) or keyword_set(nosample)) AND NOT keyword_set(zoom)) $
THEN message,"Keywords (NO)SAMPLE are effective only if ZOOM is set" 

zmx=1.

;preparation du tableau a afficher dans le cas complexe:
if datacomplex then begin
zmx=2.
if keyword_set(reim) then begin
if (reim eq 1) then arraybis = [float(array),imaginary(array)]
if (reim eq 2) then arraybis = [bytscl(float(array)),bytscl(imaginary(array))]
end else arraybis = [bytscl(abs(array)),bytscl(atan(imaginary(array),float(array)))]
end else arraybis = array

sz=size(arraybis)
nc=sz[1]
nl=sz[2]

;gestion du zoom:
if keyword_set(zoom) then begin
    
    zm=fltarr(2)
    dezoom = 0 ;si on magnifie , 1 si on de-magnifie
    
    if (zoom[0] ne -1) then BEGIN
        IF keyword_set(size) THEN $
            message,"Keyword SIZE is effective only if ZOOM=-1",/info
        zm[*]=zoom
        IF min(zm) LT 1. THEN dezoom = 1
    END else BEGIN
        IF NOT keyword_set(size) THEN size = 256
        petitcote = min([nc/zmx, nl])
        zm[*] = size/float(petitcote)
        IF petitcote GT size THEN dezoom = 1
    end
    
    if keyword_set(nosample) OR (dezoom EQ 1) then BEGIN
        ;IF NOT keyword_set(nosample) THEN message, 'Interpolation quand on dé-zoome', /info
        ncprime = fix( zm[0] * nc)
        nlprime = fix( zm[1] * nl)
        arraybis = congrid(arraybis, ncprime, nlprime, /INTERP)
    ENDIF ELSE BEGIN
        zm =  fix(zm)
        arraybis = rebin(arraybis,zm[0]*nc,zm[1]*nl,/SAMPLE)
    ENDELSE 
ENDIF ELSE zm = [1, 1]

sz=size(arraybis)

;fenetre d'affichage:
current_window = !D.WINDOW ; memorisation fenetre courante
if not keyword_set(nowin) then begin
    if keyword_set(free) then BEGIN
        win = freewnumber()
    ENDIF ELSE begin
    if NOT keyword_set(win) THEN IF current_window NE -1 THEN win = $
        current_window ELSE win = 0
    ENDELSE
    command = 'window,win,xsize=sz[1],ysize=sz[2]'
   
    IF zm[0] EQ zm[1] THEN strzm = nbr2str(zm[0]) ELSE strzm =  nbr2str(zm[0]) $
        + ' x ' + nbr2str(zm[0])  
    IF NOT keyword_set(title) THEN title = 'IDL'
;    if keyword_set(title) THEN  $
    command = command + ',TITLE='''+nbr2str(win)+': '+title+ ' ' + $
        '('+nbr2str(nc)+'x'+nbr2str(nl)+', z='+strzm+')'+''''
	if (n_elements(xpos) NE 0)  then command = command + ',XPOS=xpos'
	if (n_elements(ypos) NE 0)  then command = command + ',YPOS=ypos'
	res = execute(command)
ENDIF 

;affichage final:
if (n_elements(posinwin) ne 0) then begin
	if (size(posinwin))[0] eq 0 then $
	tvscl,arraybis,posinwin $
	else $
	tvscl,arraybis,posinwin[0],posinwin[1]
end else tvscl,arraybis

IF keyword_set(profile) THEN profiles, arraybis
IF keyword_set(rdpix) THEN rdpix,arraybis 

;dans l'option /free on repointe sur la fenetre courante:
if keyword_set(free) THEN IF (current_window GE 0) THEN wset, current_window

end







