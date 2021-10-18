; $Id: eclat.pro,v 1.2 2010-10-25 16:56:49 mvannier Exp $
function eclat,imag,INVERSE=inverse, HELP = help, VERSION = version

;+
;FONCTION
;	ECLAT - 'eclate un tableau aux quatres coins'.
;
;CATEGORIE : Array and Image Processing Routines
;
;SYNTAXE :
;  ECLAT(array)
;
;DESCRIPTION
;	ECLAT retourne l "image" eclatee au "4 coins"
;  eclat(tab) = shift(tab,-N/2,...)
;  eclat(tab,/INV) = shift(tab,+N/2,...)
;
; Remarque: justification de shift -N/2 dans le sens direct
;     N impair:
;     prenons une image 5x5 "centree" cad centre en [2,2]
;     eclater => shift de -2 = partie entiere de -N/2
;     N pair:
;     si N est pair -N/2 ou +N/2 c'est du kif!!
;  d'ou la convention adoptee.
;
;ARGUMENTS
;	ARRAY - Tableau a une,deux ou trois dimensions
;
; ATTENTION: les tableaux 3-D sont supposes etre une collection
;   de Nb images N x N  -> N x N x Nb . ECLAT eclate donc chacune de
;   ces images: shift(+-N/2,+-N/2,0) et pas shift(+-N/2,+-N/2,+-Nb/2) !!!!!!
;
;AUTEUR :
;   $Author: mvannier $
;
;
;HISTORIQUE :
;   $Log: not supported by cvs2svn $
;   Revision 1.2  2002-10-25 11:26:32+02  conan
;   les parentheses illicites ont ete remplacees par des []
;
;   Revision 1.1  2002-10-25 11:18:39+02  conan
;   Initial revision
;
;
;PREHISTORIQUE:
;	07 Juil 1993 - JM CONAN        
;  28 Juin 1994 - JM CONAN generalisation au cas 1-d et 3-D
;  05 Dec. 1995 - JM CONAN chgt des shift(tab,dim/2,..) en shift(tab,-dim/2,..)
;                 qui est plus logique pour les dim impairs
;						attention: eclat(eclat(tab)) = tab uniqt si dim pair.
;  11 Dec. 1995 - JM CONAN ajout de l'option INVERSE:
;                   eclat(eclat(tab),/INVERSE) = eclat(tab) que dim soit pair
;                   ou impair.
;-

on_error,2
IF keyword_set(version) THEN $
    printf, -2, '% '+routine_courante()+': $Revision: 1.2 $, $Date: 2010-10-25 16:56:49 $'

IF (n_params() NE 1) OR keyword_set(help) THEN BEGIN 
    message, 'Aide demandée ou syntaxe incorrecte. Documentation :', /INFO
    doc_library, routine_courante()
    retall
ENDIF

if keyword_set(inverse) then sens = 1 else sens = -1
sszz= size(imag)

nc=sszz[1]
nl=sszz[2]

if (sszz[0] ne 1) and (sszz[0] ne 2) and (sszz[0] ne 3) then $
message,"ECLAT NE TRAITE QUE DES TABLEAUX 1-D 2-D ou 3-D"
if sszz[0] eq 1 then gami=shift(imag,sens*sszz[1]/2) 							;tableau mono-dim
if sszz[0] eq 2 then gami=shift(imag,sens*sszz[1]/2,sens*sszz[2]/2) 		;image
if sszz[0] eq 3 then gami=shift(imag,sens*sszz[1]/2,sens*sszz[2]/2,0)	;cube d'image

return,gami
end

