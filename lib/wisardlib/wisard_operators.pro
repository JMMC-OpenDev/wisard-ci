; $Id: wisard_operators.pro,v 1.9 2010-10-25 16:56:50 mvannier Exp $ 
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
;   edited by Jérôme Idier, ISTE, London, 2008.
;

;+
; NAME:
;	(ROUTINE_NAME) - (Summary of purpose)
;
; CATEGORY:
;	Put a category here as defined in IDL doc.
;   See http://boa/idl/aq_idl.html for a list of categories
;	For example:
;	Array and Image Processing Routines 
;
; CALLING SEQUENCE:
;	ops=WISARD_OPERATORS( n_tels [,VERSION = version][, HELP = help])
;
; PURPOSE:
;   This function computes and returns a structure containing all the
;   operators, dimensions or algebraic structures needed in wisard. The
;   structure has the following fields:
;    N_TELS          INT       INT                          number of telescopes
;    N_BASES         INT       INT                          number of instantaneous baselines
;    N_CLOT          INT       INT                          number of instantaneous independant
;                                                           triplets containing tel #1
;    B               DOUBLE    Array[N_BASES, N_TELS]       baseline operator
;    _B              DOUBLE    Array[N_BASES, N_TELS-1]     baseline operator
;                                                           without first row
;    T_B             DOUBLE    Array[N_TELS-1, N_BASES]     transpose(_B)
;    DAG_B           DOUBLE    Array[N_TELS-1, N_BASES]     pseudo inverse of _B
;    C               DOUBLE    Array[N_CLOT, N_BASES]       Closure operator
;    DAGC            DOUBLE    Array[N_BASES, N_CLOT]       pseudo inverse of C
;    XIC             DOUBLE    Array[N_BASES, N_CLOT]       another right
;                                                           inverse of C
;
;
;
; POSITIONAL PARAMETERS:
;	n_tels   : (input) number of telescopes
;	
; KEYWORD PARAMETERS:
;
;   /VERSION : (optional input) prints version number before execution.
;   
;   /HELP    : (optional input) prints the documentation and exits.
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
;   edited by Jérôme Idier, ISTE, London, 2008.
;
; EXAMPLE:
;	ops = WISARD_OPERATORS(6)
;
; SEE ALSO:
;   All the WISARD_*.pro files, which are part of WISARD and covered by the
;   same copyright and license.
;   
; HISTORY:
;   $Log: not supported by cvs2svn $
;   Revision 1.7  2008/01/08 16:40:52  mugnier
;   Fixed typo.
;
;   Revision 1.6  2007/10/31 12:44:55  mugnier
;   Added 3rd reference.
;
;   Revision 1.5  2007/01/29 10:28:43  meimon
;   doc added
;
;   Revision 1.4  2007/01/09 11:03:55  meimon
;   remplacement de dagB par dag_B qui lui est defini (evite arithetic error,
;   devided by 0)
;
;   Revision 1.3  2006/11/02 11:03:01  meimon
;   ajout de t_B
;
;   Revision 1.2  2006/10/31 18:30:37  meimon
;   casté en double, test sur BIC 04' ok
;
;   Revision 1.1  2006/10/30 12:42:44  meimon
;   Initial revision
;
;-






; $Id: wisard_operators.pro,v 1.9 2010-10-25 16:56:50 mvannier Exp $ 
FUNCTION WISARD_OPERATORS, n_tels, VERSION = version, HELP = help

on_error,2
IF keyword_set(version) THEN $
    printf, -2, '% '+routine_courante()+': $Revision: 1.9 $, $Date: 2010-10-25 16:56:50 $'

IF (n_params() NE 1) OR keyword_set(help) THEN BEGIN ;(fill-in from syntax)
    message, 'Help required or incorrect syntax. Documentation:', /INFO
    doc_library, routine_courante()
    retall
ENDIF


n_bases = n_tels*(n_tels-1)/2
n_clot = (n_tels-1)*(n_tels-2)/2
B =  WISARD_MAT_B(n_tels)
D = WISARD_MAT_D(n_tels)
temp = WISARD_MAT_B(n_tels-1)
_B = B[*, 1:*]
C = [[-temp],[identity(n_clot)]]
iCtC = invert(matrix_multiply(C,C,/BTRANSPOSE))
ind = where(abs(iCtC) LT 1e-5, count)
IF count NE 0 THEN iCtC[ind] = 0
dagC = matrix_multiply(C,iCtC,/ATRANSPOSE) 
dag_B = matrix_multiply(invert(matrix_multiply(_B,_B,/ATRANSPOSE)),_B,/BTRANSPOSE)
xiC = transpose([[-temp*0],[identity(n_clot)]])
operators = {n_tels:n_tels, n_bases:n_bases, n_clot:n_clot, $
             B:B*1D, _B:_B*1D, t_B:transpose(_B*1D), dag_B:dag_B*1D, $
             C:C*1D, dagC:dagC*1D, xiC:xiC*1D, D:D*1D}
return, operators
END
ops = WISARD_OPERATORS(6)

END
