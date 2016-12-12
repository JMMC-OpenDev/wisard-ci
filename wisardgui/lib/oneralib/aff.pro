; $Id: aff.pro,v 1.2 2010-10-25 16:56:49 mvannier Exp $ ;(Mettre en 1ere ligne)
PRO aff,imag,pos,SAMPLE=sample,NOSAMPLE=nosample,ZOOM=zoom,SIZE=size,REIM=reim, $
    PROFILE = profile, RDPIX = rdpix, HELP = help, VERSION = version

;+
;PROCEDURE
;	AFF - version amelioree de TVSCL avec zoom...
;
;CATEGORIE
;  Image Display Routines
;
;SYNTAXE
;  AFF,array,pos,SAMPLE=sample,NOSAMPLE=nosample,ZOOM=zoom,REIM=reim
;
;DESCRIPTION
;	AFF routine d affichage avec zoom        
;
;ARGUMENTS
;	ARRAY - Tableau a deux dimensions
;	pos    - position d affichage (no de position ou [posx,posy])
;  /SAMPLE - pour zoom par macro-pixels (obsolete car = defaut)
;/NOSAMPLE - pour zoom par interpolation (devient par defaut qd on
;            de-magnifie)
;	ZOOM  - zoom: soit un entier => zx=zy=zoom ou intarr(2) => [zx,zy]=zoom
;               par defaut zoom automatique au plus proche de SIZExSIZE
;   SIZE  - taille pour le zoom automatique (par defaut 256) [cf. ZOOM] 
;	REIM  - affiche Re(array)/Im(array) si array est un tableau complexe
;              a la meme echelle   si REIM=1
;              en pleine dynamique si REIM=2
;          (par defaut on affiche MODULE/PHASE)
;  PROFILES  - lance profiles dans la foulee
;  RDPIX     - lance rdpix dans la foulee
;   /VERSION : (entrée) affichage de la version avant l'exécution.
;   
;   /HELP    : (entrée) affichage de la syntaxe et sortie du programme.
;
;VOIR AUSSI
; BYTSCL (idl) ,TVSCL (idl), TVWIN (jmc), ZTV (eric)
;
;AUTHOR
;   $Author: mvannier $
;
;HISTORIQUE:
;   $Log: not supported by cvs2svn $
;   Revision 1.1  2001-06-12 18:55:39+02  conan
;   Initial revision
;
;	07 Juil 1993 - JM CONAN        
;	31 Nov. 1995 - JM CONAN par defaut on valide l'option /sample
;                          pour interpoler => /nosample
;                          je laisse le mot clef sample qui ne sert
;                          plus a grand chose pour la compatibilite.
;  08 Jan. 1996 - JM CONAN aff n'est plus qu'un appel a la fonction
;                          plus puissante tvwin.
;  30 Avr. 1996 - JMC      ajout option REIM (pour affichage complexes)
;  13 Jun. 1998 - JMC      ajout option SIZE
;-

on_error,2
IF keyword_set(version) THEN $
    printf, -2, '% '+routine_courante()+': $Revision: 1.2 $, $Date: 2010-10-25 16:56:49 $'

IF  ((n_params() ne 1) and (n_params() ne 2)) OR keyword_set(help) THEN BEGIN 
    message, 'Aide demandée ou syntaxe incorrecte. Documentation :', /INFO
    doc_library, routine_courante()
    retall
ENDIF
 
sz = size(pos)
if not ( (sz[0] eq 0) or ((sz[0] eq 1 ) and (sz[1] eq 2 )) ) then $
message,"pos= no ou pos= [X,Y]"

if not keyword_set(imag) then message,"argument must be an 2-D array"

if not keyword_set(nosample) then nosample = 0
;dans aff.pro comme dans tvwin.pro l'option macro-pixels [sample] est prise par defaut

if not keyword_set(zoom) then zoom = -1
;par convention zoom=-1 veut dire au plus proche de 256x256 dans tvwin

if (n_params() eq 1) then $
	tvwin,imag,NOSAMPLE=nosample,ZOOM=zoom,SIZE=size,/nowin,REIM=reim, PROFILE = profile, RDPIX = rdpix $
else $ $
    tvwin,imag,NOSAMPLE=nosample,ZOOM=zoom,SIZE=size,/nowin,POSINWIN=pos,REIM=reim, PROFILE = profile, RDPIX = rdpix

end
