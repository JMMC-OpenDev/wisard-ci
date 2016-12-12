; $Id: wisard_data2mdata.pro,v 2.6 2010-10-25 16:56:50 mvannier Exp $ 
; This file is part of the WISARD software.
;
; Copyright (c) S. Meimon and L. Mugnier, ONERA, 2004-2007.
;
; This software, WISARD, is made of wisard.pro and all the wisard_*.pro files
; under the same directory, as well as the example batch file wisard_batch.pro
; 
; This software is copyright ONERA.
; The authors are Serge Meimon and Laurent Mugnier.
; E-mail: meimon at onera.fr and mugnier at onera.fr
;
; This software is a computer program whose purpose is to reconstruct an image
; of an observed object from a set of interferometric measurements.
;
; This software is governed by the CeCILL-B license under French law and
; abiding by the rules of distribution of free software. You can use, modify
; and/ or redistribute the software under the terms of the CeCILL-B license as
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
; knowledge of the CeCILL-B license and that you accept its terms.
;
; Additionally, if you use this code or a code derived from it (e.g., for a
; publication), please cite the following papers, on which it is based:
;   
; (1) S. Meimon, L. M. Mugnier, and G. Le Besnerais,
;   "Reconstruction method for weak-phase optical interferometry", 
;   Opt. Lett., 30(14):1809-1811, July 2005.
;   PDF is on-line at: http://laurent.mugnier.free.fr/publis/Meimon-OL-05.pdf
; (2) S. Meimon, L. M. Mugnier, and Guy Le Besnerais,
;   "A convex approximation of the likelihood in optical interferometry",
;   J. Opt. Soc. Am. A, November 2005.
;   PDF is on-line at: http://laurent.mugnier.free.fr/publis/Meimon-JOSAA-05.pdf
; (3) L. M. Mugnier, G. Le Besnerais, and S. Meimon, 
;   "Inversion in optical imaging through atmospheric turbulence", 
;   chapter 10 of Bayesian Approach for Inverse Problems, 
;   edited by J�r�me Idier, ISTE, London, 2008.;
;

;+
; NAME:
;	WISARD_DATA2MDATA - Computes myopic data from input data
;
; CATEGORY:
;	Array and Image Processing Routines 
;
; CALLING SEQUENCE:
;	mdata = WISARD_DATA2MDATA( data, OPERATORS = operators[, GUESS=guess,MATH=math][, /VERSION] [, /HELP])
;
; PURPOSE:
;   This function computes and returns myopic data from input
;   data. The myopic data are aranged in a array of Structures
;   containing the following fields (see doc scientifique):
;	 VISAMP          DOUBLE    Array[N_BASES]
;    VISAMPERR       DOUBLE    Array[N_BASES]
;    VISPHI          DOUBLE    Array[N_BASES]
;    VISPHIERR       DOUBLE    Array[N_BASES]
;    FREQS_U         DOUBLE    Array[N_BASES]
;    FREQS_V         DOUBLE    Array[N_BASES]
;
;
; POSITIONAL PARAMETERS:
;	data: (input)Array of Structures containing the following fields (see doc scientifique): 
;    VIS2            DOUBLE    Array[N_BASES]
;    VIS2ERR         DOUBLE    Array[N_BASES]
;    CLOT            DOUBLE    Array[N_CLOTS]
;    CLOTERR         DOUBLE    Array[N_CLOTS]
;    FREQS_U         DOUBLE    Array[N_TELS-1]
;    FREQS_V         DOUBLE    Array[N_TELS-1]
;	
; KEYWORD PARAMETERS:
;
;   MATH     : (optional input) only if 3T. Together with guess, launches
;              wisard_3t_find_phi.pro wich finds the best starting alpha
;              compatible with the guess.
;   GUESS    : (optional input)
;   /VERSION : (optional input) prints version number before execution.
;   
;   /HELP    : (optional input) prints the documentation and exits.
;   /SIMULATED_DATA : (optional input) Use a different formula to
;   convert the data to myopic visibilities, seemingly more adapted to
;   simulations where the closures and v2 are computed with maths, not
;   with real telescopes! 
;
;
; AUTHORS:
; Serge Meimon and Laurent Mugnier (ONERA). 
; Contact: lastname at onera.fr
;
;
; RESTRICTIONS:
;   This code is copyright (c) S. Meimon and L. Mugnier, ONERA, 2004-2007. 
;   
;   WISARD is governed by the CeCILL-B license under French law and abiding by
;   the rules of distribution of free software. You can use, modify and/ or
;   redistribute the software under the terms of the CeCILL-B license as
;   circulated by CEA, CNRS and INRIA at the following URL:
;   "http://www.cecill.info".
;   See source code for the full notice.
;
;   Additionally, if you use this code or a code derived from it (e.g., for a
;   publication), please cite the following papers, on which it is based:
;   
;   (1) S. Meimon, L. M. Mugnier, and G. Le Besnerais,
;   "Reconstruction method for weak-phase optical interferometry", 
;   Opt. Lett., 30(14):1809-1811, July 2005.
;   PDF is on-line at: http://laurent.mugnier.free.fr/publis/Meimon-OL-05.pdf
;   (2) S. Meimon, L. M. Mugnier, and Guy Le Besnerais,
;   "A convex approximation of the likelihood in optical interferometry",
;   J. Opt. Soc. Am. A, November 2005.
;   PDF is on-line at: http://laurent.mugnier.free.fr/publis/Meimon-JOSAA-05.pdf
;   (3) L. M. Mugnier, G. Le Besnerais, and S. Meimon, 
;   "Inversion in optical imaging through atmospheric turbulence", 
;   chapter 10 of Bayesian Approach for Inverse Problems, 
;   edited by J�r�me Idier, ISTE, London, 2008.;
;
; EXAMPLE:
;	mdata = WISARD_DATA2MDATA(data, OPERATORS = operators, VERSION = version)
;
; SEE ALSO:
;   All the WISARD_*.pro files, which are part of WISARD and covered by the
;   same copyright and license.
;   
; HISTORY:
; $Log: not supported by cvs2svn $
; Revision 2.3  2008/01/08 11:49:03  mugnier
; Fixed typo.
;
; Revision 2.2  2007/10/31 12:40:40  mugnier
; Added 3rd reference.
;
; Revision 2.1  2007/10/30 10:20:06  meimon
; includes the 3T case
;
; Revision 2.0  2007/10/01 14:23:55  meimon
; *** empty log message ***
;
; Revision 1.7  2007/10/01 14:20:18  meimon
; includes 3T exhaustive search of a best alpha init
;
; Revision 1.6  2007/02/02 16:35:39  meimon
; translation of french comments in english
;
; Revision 1.5  2007/01/29 13:07:53  meimon
; doc added
;
; Revision 1.4  2007/01/25 14:11:01  meimon
; 5.92->6
; 1.71..->sqrt(3)
;
; Revision 1.3  2007/01/23 16:45:38  mugnier
; - Si     data.vis2 < 0                 on prend mdata.visamp=0
; - Si 0 < data.vis2 < data.vis2err/5.92 on prend mdata.visamp= sqrt(data.vis2err/2.)/1.721
; - Si data.vis2 < data.vis2err          on prend mdata.visamperr=sqrt(data.vis2err)/2
;
;-



FUNCTION WISARD_DATA2MDATA, data, OPERATORS = operators, MATH = H, GUESS = guess, SIMULATED_DATA=SIMULATED_DATA, VERSION = version, HELP = help
on_error,2
IF keyword_set(version) THEN $
    printf, -2, '% '+routine_courante()+': $Revision: 2.6 $, $Date: 2010-10-25 16:56:50 $'

IF (n_params() NE 1) OR keyword_set(help) THEN BEGIN 
    message, 'Help required or incorrect syntax. Documentation:', /INFO
    doc_library, routine_courante()
    retall
ENDIF


n_tels = n_elements(data[0].freqs_u)+1
operators = WISARD_OPERATORS(n_tels, VERSION = version)
n_dat = n_elements(data)
n_bases = operators.n_bases
n_clot=(n_tels-1)*(n_tels-2)/2 ; not n_bases/3!

mdata = replicate({visamp:dblarr(n_bases), visamperr:dblarr(n_bases), $
                   visphi:dblarr(n_bases), visphierr:dblarr(n_bases), $
                   freqs_u:dblarr(n_bases), freqs_v:dblarr(n_bases)}, n_dat)

;;; Computation of spatial frequencies
mdata.freqs_u = operators._B#data.freqs_u
mdata.freqs_v = operators._B#data.freqs_v
if (n_tags(data) eq 14) then begin ; return with visdata and viserr instead of computing blind values!
   mdata.visamp=abs(data.visdata)
   mdata.visphi=atan(data.visdata,/phase)
   mdata.visamperr=abs(data.viserr)
   mdata.visphierr=atan(data.viserr,/phase)
   return,  mdata
endif

 mdata.visamp=sqrt(0.5*(data.vis2+(data.vis2^2+2*data.vis2err^2)^0.5))
 mdata.visamperr=1.0/sqrt(1.0/(mdata.visamp^2)+2*(3*mdata.visamp^2-data.vis2)/data.vis2err^2)
  
  
;;;Compute visphi
mdata.visphi = operators.dagC#data.clot

;;;cas 3T: looks for proper integer ambiguities ;;
IF (keyword_set(H) AND keyword_set(guess) AND (operators.n_tels EQ 3)) THEN BEGIN
   print, '3T dans data2mdata'
   x = reform(guess, n_elements(guess))
   achix = reform(H#x, operators.n_bases, n_elements(mdata))
   
   struct = replicate({resphi:dblarr(n_bases), resphierr:dblarr(n_bases), $
                       achix:dcomplexarr(n_bases), visy:dcomplexarr(n_bases)}, n_dat)
   arg_hx = angle(achix)
   dataphi = mdata.visphi;operators.dagC#data.clot
   struct.resphi = arg_hx-dataphi
   struct.resphierr = mdata.visphierr
   struct.achix = achix
   struct.visy = mdata.visamp*exp(!DI*dataphi)


   opt_phi = wisard_3t_find_phi(struct, operators)

   mdata.visphi = operators.dagC#data.clot+operators._B#opt_phi.phi
ENDIF

; simulated data (beauty contests, Aspro...). Use historical formula below.
if keyword_set(simulated_data) then begin 
   mdata.visphierr = 3.*((operators.dagC)^2#data.cloterr) 
endif else begin
; else stick to interferometric observables obtention, use the  
; following formula in Tatulli et Chelli, 2005, JOSA Vol 22 no 8 Annex 1
; with corrections of multiple typos...
; we thought to estimate visphierr by reversing the estimate of sigma2_clot
; as a function of sigma2_visphi and vis2 
; possibly factor 2 (0.5 below) wrong
; use visamp = sqrt(v2) but never negative:
   denom=operators.d#(0.5/(mdata.visamp)^2)
   weightedCloterr2=data.cloterr^2/denom
   visphierr2 = 3.*((operators.dagC)^2#weightedCloterr2)/((mdata.visamp)^2)
   mdata.visphierr = sqrt(visphierr2)
endelse
return,  mdata
END

