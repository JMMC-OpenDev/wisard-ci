; $Id: wisard.pro,v 3.4 2009/11/10 $ 
; This file is part of the WISARD software.
;
; Copyright (c) S. Meimon and L. Mugnier, ONERA, 2004-2007.
; Latest adds and corrections by M. Vannier and A.Domiciano, Fizeau/UNS/OCA, 2010-2011.

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
;   edited by J�r�me Idier, ISTE, London, 2008.
;
  

;+
; NAME:
;	WISARD - Image reconstruction from interferometric data
;
; CATEGORY:
;	Array and Image Processing Routines 
;
; CALLING SEQUENCE:
;
;   IMAGE = WISARD(data,     $
;                            FOV = fov, NP_MIN = np_min, $
;                            OVERSAMPLING = oversampling, $
;                            GUESS = guess, $
;                            NBITER = nbiter, THRESHOLD = threshold, $
;                            POSITIVITY = positivity, $
;                            PSD = psd, MEAN_O = mean_o,
;                            DELTA = delta, SCALE = scale, WHITE =
;                            white, TOTVAR=TOTVAR, $
;                            MU_SUPPORT = mu_support, FWHM = fwhm, $
;                            LIBRARY = library, $ $
;                            AUX_OUTPUT = aux_output, DISPLAY = display, $
;                            [,/CHI2CRIT ] [, /VERSION] [, /HELP])
;
; PURPOSE:
; 
;   The WISARD function computes and returns the IMAGE reconstructed from the
;   interferometric input data. WISARD stands for ``Weak-phase Interferometric
;   Sample Alternating Reconstruction Device''.
;
;   It is based on (and well-described in) the PhD thesis work of S. Meimon,
;   available on the French national multidisciplinary thesis server TEL at
;   http://tel.archives-ouvertes.fr/docs/00/05/44/98/PDF/these_finale.pdf
;
;   This routine and all the WISARD_* routines necessary to it are a complete
;   rewrite (from scratch) by S. Meimon and L. Mugnier of earlier programs of
;   the aforementioned thesis.
;
;   The reconstruction is iterative and should be stopped only when the
;   reconstructed object (or the minimized criterion) no longer evolves to
;   some given precision; the latter is given through a THRESHOLD detailed
;   below.
;   However, for "quick-look" results, there are two other ways to stop the
;   iterations: 
;   - give a maximum number of iterations NBITER (see below);
;   - interrupt manually the iterations, by typing 'Q' (capital q character)
;     in the IDL interpreter.
;
;
; POSITIONAL PARAMETERS:
;
;   data         : (input) interferometric input data for
;                  the reconstruction. The data format is described in the
;                  documentation that accompanies WISARD (doc/ directory).
;                  In short, "data" is a vector of structures, one structure
;                  per time of measurement, each containing the following
;                  fields: 
;                  VIS2 : the squared visibilities for each baseline
;                  VIS2ERR : the standard deviation on the squared
;                            visibilities 
;                  CLOT : the closure phases for each triplet of telescopes
;                         involving telescope 1
;                  CLOTERR : the standard deviation on the closure phases
;                  FREQS_U : first coordinate of the spatial frequencies
;                            involving telescope 1 
;                  FREQS_V : second coordinate of the spatial frequencies
;                            involving telescope 1.
;
;                  The zero frequency must not be present in the data: the
;                  data must be normalized, as are OIFITS files, so that the
;                  value of the data at this frequency would be 1. In other
;                  words, the reconstructed object has a unit sum.           
;         
; KEYWORD PARAMETERS:
;   
;   FOV          : (input) Field-Of-View of the reconstructed image, in units
;                  consistent with the data. More precisely the unit for FOV
;                  must be the inverse of the unit of the arrays of frequencies
;                  (FREQS_U and FREQS_V) of the data. Usually FREQS_U and
;                  FREQS_V are in rd^(-1) so FOV should be in rd (radians).
;                  
;   NP_MIN       : (input) MINimum width (Number of Points) of the
;                  reconstructed image. The Number of Points of the
;                  reconstructed object may be greater, depending on FOV,
;                  OVERSAMPLING factor, and frequencies present in the data.
;                  See routine WISARD_MAKE_H for details. 
;                  
;   OVERSAMPLING : (optional input) oversampling factor for the reconstructed
;                  image. By default, OVERSAMPLING=1 and the maximum spatial
;                  frequency of the reconstructed object is the maximum
;                  frequency of the data. See routine WISARD_MAKE_H for details.
;                  
;   GUESS        : (optional input) initial guess for the reconstructed image
;                  (or dirty map if not present). This guess is massaged in
;                  the following way: resampled if necessary to the correct
;                  size, thresholded to positive values and normalized to a
;                  unit sum. This massaged guess is available on exit as a
;                  field of AUX_OUTPUT.
;                  
;   NBITER       : (optional input) maximum NumBer of ITERations for the
;                  reconstruction, 500 by default. For a better control of
;                  the reconstruction, one should rather use THRESHOLD below
;                  and not lower this value.
;                  
;   THRESHOLD    : (optional input) convergence THRESHOLD to be used as a
;                  stopping criterion for the iterations. By default set to
;                  the machine precision in simple precision ~1.19e-07
;                  (although computations are done in double precision).
;                  For a (rather) quick-look result, set to a smaller
;                  value, e.g., 1e-6. 
;                  
;   POSITIVITY   : (optional input) POSITIVITY constraint for the
;                  reconstruction. It is set to true (1) by default. Set it
;                  explicitly to 0 if by misplaced curiosity you do not want
;                  to use the positivity constraint in the reconstruction.
;                  
;   PSD          : (optional input) 2-D map of size NP_MIN x NP_MIN containing
;                  the PSD for the (quadratic) regularization of the
;                  reconstruction. See routine J_PRIOR_GAUSS and example file
;                  for details. 
;                  
;   MEAN_O       : (optional input) 2-D map of size NP_MIN x NP_MIN containing
;                  the MEAN Object to be used for regularization of the
;                  reconstruction. MEAN_O is used both for PSD regularization
;                  and for white L1-L2 regularization (i.e., when DELTA, SCALE
;                  and WHITE are set). It is also used for soft support
;                  regularization if FWHM is not given; in this case it must
;                  be a 2D map with stricly positive values.
;                  
;   DELTA        : (optional input) scalar factor for L1-L2 regularization,
;                  used to set the threshold between quadratic (L2) and linear
;                  (L1) regularization.
;                  See J_PRIOR_L1L2 and the example file for some more details.
;                  See the following paper for a complete description of the
;                  L1-L2 regularization: L. M. Mugnier, T. Fusco, and J.-M.
;                  Conan. "MISTRAL: a myopic edge-preserving image restoration
;                  method, with application to astronomical
;                  adaptive-optics-corrected long-exposure images", J. Opt.
;                  Soc. Am. A, 21(10):1841-1854, October 2004. The PDF of this
;                  paper is on-line at:
;                  http://laurent.mugnier.free.fr/publis/Mugnier-JOSAA-04.pdf
;                  
;   SCALE        : (optional input) scalar SCALE factor for L1-L2 regularization.
;                  should be of the order of the average object value if
;                  /WHITE is used, and of the order of the RMS object's
;                  gradient value if WHITE is not set. See J_PRIOR_L1L2 and
;                  example file for some more details. See the paper cited
;                  above for the whole story.
;                  
;   /WHITE       : (optional input) flag to switch between edge-preserving
;                  regularization (WHITE=0, default) and spike-preserving
;                  regularization (WHITE=1). In the latter case the
;                  regularization is performed independently on each pixel
;                  value, hence the flag name.
;                  
;   MU_SUPPORT   : (optional input) global factor for regularization by "soft
;                  support". See J_PRIOR_SUPPORT for more details. 
;                  See the following paper for a complete description of the
;                  soft support regularization: 
;                  Imaging with long-baseline optical interferometry, 
;                  G. Le Besnerais, S. Lacour, L. M. Mugnier, E. Thiebaut, G.
;                  Perrin and S. Meimon, IEEE Journal of Selected Topics in
;                  Signal Processing, submitted (2008). 
;                  
;   FWHM         : (optional input) FWHM for regularization by soft support.
;                  If given, a Lorentzian is computed as prior object and MEAN_O
;                  is discarded even if present. If not given, MEAN_O must be
;                  given and is used as the prior object.
;                  
;   LIBRARY      : (optional input) full path of the OptimPack library (see
;                  FMIN_OP for details), if necessary. Under Unix systems the
;                  OptimPack library  should be found automatically.
;                  
;   AUX_OUTPUT   : (optional output) structure containing various optional
;                  AUXiliary outputs, for diagnostic purposes.
;                  Its fields are:
;                    CMDATA          STRUCT    -> Convexified myopic data
;                                              structure, see scientific doc
;                    PRIOR           STRUCT    -> see WISARD_SET_REGUL
;                    CRIT_ARRAY      DOUBLE    Criterion components (data
;                                              fidelity and regularization
;                                              at convergence).
;                    X               DOUBLE    Final solution obtained (object)
;                    ALPHA           DOUBLE    Final solution obtained
;                                              (aberrations) 
;                    OPERATORS       STRUCT    -> structure, output of function
;                                              WISARD_OPERATORS 
;                    WEIGHTS_CONSTANT
;                                    STRUCT    -> see below 
;                    WEIGHTS_ALPHA   STRUCT    -> not documented 
;                    WEIGHTS_X       STRUCT    -> see below  
;                    RAD_FREQS       DOUBLE    Radial frequencies
;                    GUESS           DOUBLE    Initial guess used.
;                  
;       AUX_OUTPUT.WEIGHTS_X : Structure, field of the AUX_OUTPUT structure
;                  Its fields are:
;                    ABS_HX          array    -> visibility moduli of the
;                    reconstructed object; useful for plots, to be compared to
;                    the abs((aux_output.cmdata).vis), the visibility modulus
;                    pseudo data.
;                    ARG_HX,WX1,WX2  array    -> not documented
;
;       AUX_OUTPUT.WEIGHTS_CONSTANT : Structure,field of the AUX_OUTPUT
;                                     structure. Its fields are:
;                    RE_Y_DATA:                  real(cmdata.vis)
;                    IM_Y_DATA:                  imaginary(cmdata.vis)
;                    ABS_Y_DATA:                 abs(cmdata.vis)
;                    ARG_Y_DATA:                 angle(cmdata.vis)
;                    W11, W12, W22               see doc scientifique
;                    W_RAD:                      -> see cmdata.w_rad  
;                    W_TAN:                      -> see cmdata.w_tan  
;                    INDX_FLAG:                  integers (flag indices)
;                    
;   DISPLAY      : (optional input) if this keyword is set (>0) then the
;                  program displays the reconstructed object in window
;                  (display) and the fit to the visibilities in window
;                  (display+1) along the way.
;
;   LUT          : (optional input) Look-Up Table (a.k.a. color table in IDL
;                  jargon) to be used for display. Defaults to 13, i.e. rainbow.
;
;   /CHI2CRIT    : (optional input) casts the criterion in a chi2 style
;
;   /VERSION     : (optional input) prints version number before execution.
;   
;   /HELP        : (optional input) prints the on-line documentation and exits.
;
;  /COPYRIGHT    : (optional input) prints information about copyright and exits
;
; /VERBOSE_TEST : (optional input) prints intermediate results, mainly for control purposes
; 
; AUTHORS:
; Serge Meimon and Laurent Mugnier (ONERA). 
; Contact: lastname at onera.fr
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
;   edited by J�r�me Idier, ISTE, London, 2008.
;
; EXAMPLE:
;	See example batch file in the pro/ directory. 
;
; SEE ALSO:
;   All the WISARD_*.pro files, which are part of WISARD and covered by the
;   same copyright and license.
;
; ACKNOWLEDGMENTS:
;
;   The authors of WISARD hereby express their gratefulness to the following
;   good guys ;-) :
;   
;   �ric Thi�baut:
;   WISARD heavily uses OptimPack (by �ric Thi�baut, from CRAL) as its (neat)
;   minimization engine. WISARD assumes that the OptimPack library and its IDL
;   frontend (op_*.pro files) are properly installed on your machine.
;
;   Fr�d�ric Cassaing, Jean-Marc Conan and Jean-Fran�ois Sauvage:
;   WISARD uses a few routines by these ONERA scientists. These routines (and
;   some others by Laurent Mugnier) are distributed along with WISARD,
;   with these authors' permissions, in the lib/oneralib/ directory.
;
;   F K Knight, Kim Tolbert and W. Landsman:
;   WISARD uses legend.pro (part of NASA astro library, available at
;   http://idlastro.gsfc.nasa.gov/). 
;   
; HISTORY:
;   Revision 3.3, branche _flag_to_weight0:
;   * les infos de flag ne sont pas transmises pour la minimisation 
;   (à wisard_j*, à travers fmin_op), mais les flags "bad" sont traduits en weight=0
;   * Comparaison avec _v3.3 (avec transmission des flags)
;   
;   Revision 3.3 2009/11/10 mvannier
;   * Added support for VIS and CLOT flags reference: 
;   indx_goodflag[nBase, nTime*nLambda] is computed from the 'flag' matrices 
;   in the input data structure. It is now included in weights_constant structure, 
;   and transmitted to wisard_jtotal_alpha and FMIN_OP, in order to mask unwanted 
;   data during minimization process. 
;   * => TO BE CLEANED (comments, VERBOSE_TEST, etc...)
;   * => saved as wisard_v3.3.pro
;   
;   Revision 3.2  2008/09/11 13:57:17  mugnier
;   Added LUT keyword.
;
;   Revision 3.1  2008/04/17 10:36:44  mugnier
;   When DISPLAY is set, the windows used are now (DISPLAY) and (DISPLAY+1)
;   instead of 0 and 1 previously.
;
;   Revision 3.0  2008/03/28 15:38:17  mugnier
;   The "soft support" regularization is now supported by WISARD.
;
;   Revision 2.6  2008/01/08 16:40:00  mugnier
;   Fixed typo.
;
;   Revision 2.5  2008/01/08 10:56:49  mugnier
;   Call doc_library without call_procedure to see the copyright + fixed typos.
; 
;   Revision 2.4  2008/01/08 09:19:27  mugnier
;   Added JFS in acknowledgments (for read_params_vm.pro).
; 
;   Revision 2.3  2007/10/31 12:42:43  mugnier
;   Added 3rd reference.
; 
;   Revision 2.2  2007/10/30 10:18:11  meimon
;   minor changes
;   Revision 2.1 2007/10/01 13:23:13 meimon now
;   _3T_data2MDATA is redundant. We use data2MDATA instead, with the right
;   keywords
;
;   Revision 2.0  2007/10/01 13:19:15  meimon
;   3T case separately taken care of
;
;   Revision 0.13  2007/01/31 19:15:59  mugnier
;   The thresholds used for each minimization (in alpha or in x) and for the
;   global minimization are now consistent whatever the input THRESHOLD is.
;
;   Revision 0.12  2007/01/29 15:42:22  mugnier
;   Documentation enriched.
;
;   Revision 0.11  2007/01/24 17:58:10  mugnier
;   First documented version!
;
;   Revision 0.10  2007/01/22 15:59:43  meimon
;   plot_fit ->wisard_plot_fit.pro
;
;   Revision 0.9  2007/01/22 15:24:54  meimon
;   plot_fit module added
;
;   Revision 0.8  2007/01/19 17:53:07  mugnier
;   Am�liorations affichage.
;
;   Revision 0.7  2007/01/19 13:00:35  mugnier
;   On sauve le guess d'initialisation, si besoin r��chantillonn� et seuill� � 0
;   dans AUX_OUTPUT (et ce guess est la dirty map seuill�e � 0 s'il n'est pas
;   pass� en argument).
;
;   Revision 0.6  2007/01/18 17:05:27  mugnier
;   On exit, X is reform'ed into a square image.
;
;   Revision 0.5  2007/01/18 15:18:42  mugnier
;   INIT=1 no longer necessary with FMIN_OP version >= 1.9, thus removed.
;   ACTIVE_SET computed here, not in WISARD_SET_REGUL.
;   Removed remains of commented code.
;
;   Revision 0.4  2006/11/09 17:54:56  meimon
;   arret sur convergence d'un cycle � l'autre
;   joli resultat sur beauty contest 04
;
;   Revision 0.3  2006/11/07 16:25:43  meimon
;   version serge
;-

FUNCTION WISARD, data,  $
                 FOV = fov, NP_MIN = np_min, $
                 OVERSAMPLING = oversampling, $
                 GUESS = guess, VERBOSE_TEST=verbose_test, $
                 NBITER = nbiter, THRESHOLD = threshold, $
                 POSITIVITY = positivity, $
                 PSD = psd, MEAN_O = mean_o, $
                 DELTA = delta, SCALE = scale, WHITE = white, TOTVAR=TOTVAR, $
                 MU_SUPPORT = mu_support, FWHM = fwhm, $
                 LIBRARY = library, $
                 AUX_OUTPUT = aux_output, CHI2CRIT = chi2crit, $
                 DISPLAY = display, LUT = lut, $ 
                 VERSION = version, HELP = help, $
                 COPYRIGHT = copyright, PRINT_TIMES=print_times, USE_FLAGGED_DATA=use_flagged_data, DEBUG=debug,_EXTRA=ex
  
  ;; for tests: 
   t_fminop=0d;
  
  on_error,2
  IF keyword_set(version) THEN $
     printf, -2, '% WISARD: $Revision: 3.2 $, $Date: 2008/09/11 13:57:17 $'

  IF keyword_set(copyright) THEN BEGIN
     doc_library, 'wisard_license'
     message, 'if you accept these conditions you can now run WISARD.', /INFO
     retall
  ENDIF

  IF (n_params() NE 1) OR keyword_set(help) THEN BEGIN 
     message, 'Help required or incorrect syntax. Documentation:', /INFO
     doc_library, 'wisard'
     retall
  ENDIF

  message,'* Copyright (C) S. Meimon and L. Mugnier, ONERA, 2004-2007.', /info
  message,'This software, i.e., WISARD and the wisard_* routines used ' + $
          'by it, are Copyright ONERA.', /info
  message,'The authors of WISARD are Serge Meimon and Laurent Mugnier.', /info
  message,'The use of WISARD implies acceptance of the conditions described ' + $
          'in wisard_license.pro', /info

  IF NOT keyword_set(display) THEN display = 0B ; do not display by default

  IF n_elements(VERBOSE_TEST) eq 0 then VERBOSE_TEST=0 ;else VERBOSE_TEST=1

  if n_elements(use_flagged_data) eq 0 then use_flagged_data=0 ; 

  IF keyword_set(display) THEN BEGIN
     IF n_elements(lut) EQ 0 THEN $
        loadct, 13 $            ; display uses rainbow LUT by default
     ELSE $
        loadct, lut
  ENDIF
  
  if keyword_set(print_times) then t=SYSTIME(/SECONDS )
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; FABRICATION DONNEES COMPLEXES MYOPES et matrice H
  mdata = WISARD_DATA2MDATA(data, OPERATORS = operators, VERSION = version,_EXTRA=ex)
  cmdata = WISARD_MDATA2CMDATA(mdata, OPERATORS = operators, VERSION = version,_EXTRA=ex)

  matrix_phases=angle(cmdata.vis)
  clot_from_add_phases=operators.C#matrix_phases

  if keyword_set(print_times) THEN t=SYSTIME(/SECONDS)
  mult_phasors=multiply_phasors(cmdata.vis,operators.C) ; in wisardlib : 
  if keyword_set(print_times) THEN print,' Time multiply : ', SYSTIME(/SECONDS )-t
  
  if keyword_set(print_times) THEN t=SYSTIME(/SECONDS)
  cloture_from_cmdata=angle(mult_phasors)
  if keyword_set(print_times) THEN print,' Time angle : ', SYSTIME(/SECONDS )-t
 
  vec=cmdata.w_rad              ;
 
  IF keyword_set(display) THEN BEGIN
     IF n_elements(lut) EQ 0 THEN $
        loadct, 13 $            ; display uses rainbow LUT by default
     ELSE $
        loadct, lut
  ENDIF


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Flag indices 
  median_err=median(data.vis2err)
  flagvis=(data.vis2flag or (abs(data.vis2) gt 1.+3*median_err) or (data.vis2err le 1e-10) ) 
  
  if (~use_flagged_data) then $
     indx_goodflag=where(((flagvis EQ 0) AND (abs(operators.dagc)#(data.clotflag) EQ 0)), $
                         n_indx_goodflag, complement=indx_badflag) $
  else begin
     n_indx_goodflag=n_elements(flagvis) ;
     indx_goodflag=indgen(n_indx_goodflag)
  endelse

  n_indx_badflag=n_elements(flagvis)-n_indx_goodflag 
  
  if ((n_indx_badflag GT 0) AND (VERBOSE_TEST)) THEN begin
     indx_badflag_2d=array_indices(flagvis,indx_badflag) 
     print, ' abs(cmdata.vis)[indx_badflag] ', (abs(reform((cmdata.vis),n_elements(cmdata.vis))))[indx_badflag]
     print, n_indx_badflag,' bad-flagged elements' 
     indx=where((abs(operators.dagc)#(data.clotflag) NE 0),count)
     print,'   ', count,' from flagged closure.'
     indx=where((flagvis NE 0),count)
     print,'   ', count,' from flagged VIS2.'
  ENDIF
  if keyword_set(print_times) then begin
     print,' Time INIT: (c)mdata, flags : ', SYSTIME(/SECONDS )-t
     t=SYSTIME(/SECONDS)
  endif

;
  freqs_u=reform(cmdata.freqs_u, operators.n_bases*n_elements(mdata))  
  freqs_v=reform(cmdata.freqs_v, operators.n_bases*n_elements(mdata))
  H=WISARD_MAKE_H(FREQS_U=freqs_u, FREQS_V=freqs_v,$
                  OVERSAMPLING = oversampling,$
                  FOV = fov, NP_MIN = np_min,$
                  ORIGIN = origin, $ 
                  NP_OUTPUT = NP, STEP_OUTPUT = step_output, /verbose, VERSION = version)
                      
  if keyword_set(print_times) then begin 
     print,' Time INIT: H : ', SYSTIME(/SECONDS )-t
     t=SYSTIME(/SECONDS)
  endif
  
;;;;;;;;;;Rebuild in case NT=3, to avoid integer ambiguities
  IF ((operators.n_tels EQ 3) AND (keyword_set(guess))) THEN BEGIN
     print, '3T case'
     IF keyword_set(guess) THEN nguess = congrid(guess, NP, NP)  > 0
     mdata = WISARD_DATA2MDATA(data, OPERATORS = operators, MATH = H, guess = nguess, VERSION = version)
     cmdata = WISARD_MDATA2CMDATA(mdata, OPERATORS = operators, VERSION = version)   
     if keyword_set(print_times) then print,' Time INIT: (c)mdata for 3T case : ', SYSTIME(/SECONDS )-t
  ENDIF
  
  t=SYSTIME(/SECONDS)

;; weights set at 0 for all 'bad'flags   
  if (n_indx_badflag gt 0) then begin 
     print, format='(%" Setting a null weight on %d bad-flagged data")',n_indx_badflag
     
     w_tan_0=cmdata.w_tan       ; intermediate variable 
     w_rad_0=cmdata.w_rad       ; intermediaire variable 
     w_tan_0[indx_badflag]=(w_rad_0[indx_badflag]=0.)
     cmdata.w_tan=w_tan_0
     cmdata.w_rad=w_rad_0 
  endif
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; INITIALISATIONS

; Initialisation of default object: dirty map, or guess if given, with positivity
  IF keyword_set(guess) THEN BEGIN
     nguess = congrid(guess, NP, NP)  > 0 ; matrice guess redimenssionnee a NPxNP
  ENDIF  ELSE BEGIN
; inverse FT of myopic complex visibilities:
     if (verbose_test) then t0=systime(1)
     if (verbose_test) then print, "dirty map memory at start: ",memory(/current)
     dirty_map=0
     a=(conj(H))[indx_goodflag,*]
     b=(reform(cmdata.vis,(size(H))[1]))[indx_goodflag]
     dirty_map=matrix_multiply(a,b,/ATRANSPOSE)
     a=0
     b=0
     if (verbose_test) then print,'time making dirty map:', systime(1)-t0
     if (verbose_test) then print, "memory at end: ",memory(/current)
     nguess = real(dirty_map, VERSION = version) > 0
  ENDELSE
  nguess = nguess/total(nguess) ; saved in aux_output.
  x = reform(nguess, n_elements(nguess))

  alpha = randomn(seed, operators.n_tels-1, n_elements(cmdata))*0D ;

;Use pointers for a faster passage of variables
  re_y_data = PTR_NEW(real(cmdata.vis), /NO_COPY)
  im_y_data = PTR_NEW(imaginary(cmdata.vis), /NO_COPY)
  abs_y_data = abs(cmdata.vis)
  arg_y_data = angle(cmdata.vis)
  w11 = PTR_NEW(cmdata.w_rad*abs2(COS(arg_y_data)) $
                +cmdata.w_tan*abs2(SIN(arg_y_data)), /NO_COPY)
  w22 = PTR_NEW(cmdata.w_tan*abs2(COS(arg_y_data)) $
                +cmdata.w_rad*abs2(SIN(arg_y_data)), /NO_COPY)
  w12 = PTR_NEW((cmdata.w_rad-cmdata.w_tan)* $
                COS(arg_y_data)*SIN(arg_y_data), /NO_COPY)

  weights_constant = {re_y_data:re_y_data, $
                      im_y_data:im_y_data, $
                      abs_y_data:abs_y_data, $
                      arg_y_data:arg_y_data, $
                      w11:w11, $
                      w12:w12, $
                      w22:w22, $ 
                      w_rad:(cmdata.w_rad), $ 
                      w_tan:(cmdata.w_tan)}
  
  rad_freqs = reform(sqrt(abs2(cmdata.freqs_u)+abs2(cmdata.freqs_v)), $
                     operators.n_bases*n_elements(cmdata))
  clot_freqs=cmdata_clot_maxfrequency(cmdata)

  if keyword_set(debug) then begin 

;;;;;;;;;;;;;
;; plot of V2 from data and from myopic convex estimate
     window, display+3, xs = 512, ys = 512, $
             title = 'IDL '+nbr2str(display)+': original and myopic vis.'  
     wset, display+3
     RAD_FREQS = rad_freqs                                     ;
     ABS_VIS2 = reform(sqrt(abs(data.vis2)), n_elements(data.vis2)) ;
     ABS_Y = reform(weights_constant.abs_y_data, n_elements(ABS_VIS2))
     INV_SIGMA_RAD = reform(sqrt(abs(data.vis2))*0.0, n_elements(ABS_VIS2))

     xtitle =  'frequency'      ; textoidl('|\nu|') 
     ytitle = 'VIS.'
     xmax=max(rad_freqs)*1.05D
     xmin=0.0D

     !P.MULTI=[0,0,2,0,0]
     !P.POSITION=[0.05,0.3,0.95,0.95]
     plot,RAD_FREQS,ABS_VIS2, YRANGE=[0,1], XRANGE = [xmin, xmax], XMARGIN=[10,4], /NODATA ; white, box only
     oplot,rad_freqs,abs_VIS2,psym=7, color = 150 ; red
     oplot, rad_freqs, psym=5, reform(mdata.visamp,n_elements(ABS_VIS2)), color=60
     oplot,rad_freqs,abs_y,psym=6,color = 100 ;green

     tstring = ['Original sqrt(vis2)', 'Convex (mdata)','Myopic convex (cmdata)']
     tlinestyle= [0,0,0]        ;
     tsym= [7, 5, 6]
     wis_legend, tstring, line = tlinestyle, psym = tsym, color = [150,60,100], /top, $
                 /right, clear = clear
     !P.MULTI=0
     !P.POSITION=0
  endif
  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; REGULARISATION
  prior = WISARD_SET_REGUL( POSITIVITY = positivity, $
                            PSD = psd, MEAN_O = mean_o, $
                            MU_SUPPORT = mu_support, FWHM = fwhm, $
                            SCALE = scale, DELTA = delta, WHITE = white,  $
                            TOTVAR = TOTVAR, $
                            NP = NP, VERSION = version)
  IF (prior.positivity NE 0B) THEN active_set = bytarr(prior.squareNP) + 1B

  if keyword_set(print_times) then begin 
     print,' Time INIT: init x, w, w_constant, H, regul : ', SYSTIME(/SECONDS )-t
     t=SYSTIME(/SECONDS)
  endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; MAIN LOOP (SELF-CALIBRATION)
  continue = 1B
  iter = 0L                     ; idem for fmin_gpa et fmin_op
  crit_array = dblarr(2)        
  epsilon = (machar()).eps      ; simple-precision ~1.19e-07

  IF keyword_set(threshold) THEN epsilon = threshold
  IF NOT keyword_set(nbiter) THEN nbiter = 500L

  itmax_alpha = operators.n_tels ; -1L ; NTELS is needed here
  itmax_x = NP/12L+3L            ;itmax_x = prior.squareNP

  conv = epsilon+1D
  IF (display NE 0) THEN BEGIN
     window, display, xs = 256, ys = 256, $
             title = 'IDL '+nbr2str(display)+': reconstruction'
     window, display+1, $
             title = 'IDL '+nbr2str(display+1)+': plot of visibility fit'
  ENDIF

  if keyword_set(print_times) then print, '*** Starting loop'

  WHILE ((iter LT nbiter) AND (conv GE epsilon) AND (continue)) DO BEGIN
     
     if keyword_set(print_times) then t_loop= SYSTIME(/SECONDS ) 
     if (VERBOSE_TEST) THEN print, ' ********* iter : ', iter
     old_crit = total(crit_array)
     ;; MIN sur alpha
     ;;maj de weights_x
     factor = 1D/total(x)
     norm_x = x*factor
     
     if (VERBOSE_TEST) THEN print, ' rms(map x) : ',stddev(x)
     
     if keyword_set(print_times) then t= SYSTIME(/SECONDS )
     achix = reform(H#norm_x, operators.n_bases, n_elements(cmdata)) ; TF de norm_x (guess au depart, resultat de boucle ensuite)  
     if keyword_set(print_times) then print,'  Time: achix orig: ', SYSTIME(/SECONDS )-t
      
     abs_hx = abs(achix) 
     abs2_hx = abs2(achix)
     arg_hx = angle(achix)

     wx1 = abs2_hx*(cmdata.w_tan-cmdata.w_rad)
     wx2 = 2D*abs_hx*abs_y_data*cmdata.w_rad

     weights_x = {abs_hx:abs_hx, arg_hx:arg_hx, wx1:wx1, wx2:wx2}
     
     IF (display NE 0) THEN BEGIN
        if keyword_set(print_times) then t= SYSTIME(/SECONDS ) ; t4 
        wset, display+1

    mult_phasors=multiply_phasors(achix,operators.C)
    cloture_from_current_x=angle(mult_phasors)          

        WISARD_PLOT_FIT, RAD_FREQS = rad_freqs,$ 
                         ABS_HX = reform(abs_hx, n_elements(abs_hx)), $
                         ABS_Y  = reform(weights_constant.abs_y_data, n_elements(abs_hx)), $
                         CLOT_FREQS = clot_freqs,$
                         CLOT_FROM_DATA = reform(data.clot, n_elements(data.clot)),$
                         CLOT_FROM_CMDATA = reform(cloture_from_cmdata, n_elements(cloture_from_cmdata)),$
                         CLOT_FROM_CURR_X = reform(cloture_from_current_x, n_elements(cloture_from_cmdata)),$
                         INV_SIGMA_RAD = reform(sqrt(cmdata.w_rad),n_elements(abs_hx)), $
                         INDX_BADFLAG=indx_badflag,$ ; uncomment if bad flagged data is to be plotted
                         VERSION = version
        if keyword_set(print_times) then print,'  Time: Display :',SYSTIME(/SECONDS )-t
     ENDIF
     
     ;;minimization of alpha    
                                ;
     if keyword_set(print_times) then t= SYSTIME(/SECONDS ) ; t4 
     FMIN_OP, alpha, error, FUNC = 'wisard_jtotal_alpha', LIBRARY = library, $
              /ALLTHEWAY, CONV_THRESHOLD = epsilon/2D, $
              ITMAX = itmax_alpha, $
              _BALPHA = _Balpha, $
              WEIGHTS_CONSTANT = weights_constant, WEIGHTS_X = weights_x, $
              CRIT_ARRAY = crit_array, $
              OPERATORS = operators, CHI2CRIT = chi2crit, VERSION = version
     if keyword_set(print_times) then print,'  Time: FMIN_OP(alpha) : ',SYSTIME(/SECONDS )-t
     
     ;; MIN of x, possibly with positivity constraint   
     
     ; build ph matrix (fast method): 
     if keyword_set(print_times) then t= SYSTIME(/SECONDS )                                 
     ph=H*(REFORM(EXP(!dI*(operators._B#alpha)),operators.n_bases*n_elements(cmdata),/overwrite)#replicate(1.,n_elements(H[0,*])))
     if keyword_set(print_times) then print, '  Time: PH=P.H : ',SYSTIME(/SECONDS )-t
     
     ;Use pointers for a faster passage of variables
     weights_alpha = {re_ph:PTR_NEW(real(ph), /NO_COPY), im_ph:PTR_NEW(imaginary(ph), /NO_COPY)}
     
     ;; minimisation over X
     if keyword_set(print_times) then t= SYSTIME(/SECONDS ) 
     FMIN_OP, x, error, FUNC = 'wisard_jtotal_x', LIBRARY = library, $
              /ALLTHEWAY, CONV_THRESHOLD = epsilon/2D, $
              ITMAX = itmax_x,  $
              ACTIVE_SET = active_set, $
              WEIGHTS_CONSTANT = weights_constant, $
              WEIGHTS_ALPHA = weights_alpha, $
              PRIOR = prior, CRIT_ARRAY = crit_array, $
              OPERATORS = operators, CHI2CRIT = chi2crit, VERSION = version
                   
     if keyword_set(print_times) then BEGIN
      t_fminop=t_fminop+SYSTIME(/SECONDS )-t
      print, '   TIME: FMIN_OP(x) : ',SYSTIME(/SECONDS )-t
      t= SYSTIME(/SECONDS ) 
     endif 
     
     IF (display NE 0) THEN BEGIN
        wset, display                 
        aff, rotate(reform(x, NP, NP), 5), VERSION = version
     ENDIF
     if keyword_set(print_times) then print,'Time t8',SYSTIME(/SECONDS )-t
     
     iter = iter+1L
     IF (getenv('DISPLAY') NE '') THEN BEGIN 
        IF (strupcase(get_kbrd(0) EQ 'Q')) THEN BEGIN   ; Q: clean interruption of iterations
           continue = 0B
           message, 'interruption of iterations by user.', /INFO
        ENDIF
     ENDIF
     IF iter GT 2 THEN BEGIN
        conv = 2D*(old_crit-total(crit_array))/(old_crit+total(crit_array))
        print,nbr2str(iter, VERSION = version), '/', nbr2str(nbiter), '. Convergence=',  $
              nbr2str(conv, format = '(E25.2)'), $
              '. Criterion=',  $
              nbr2str(total(crit_array), format = '(F25.4)'), ' = ', $
              nbr2str(total(crit_array[0]), format = '(F25.4)'), ' + ', $
;              nbr2str(total(crit_array[0])/(2*n_elements(cmdata)), format = '(F25.4)'), ' + ', $
              nbr2str(total(crit_array[1]), format = '(F25.4)')
     ENDIF
     if keyword_set(print_times) then print,'Total time for iteration : ',SYSTIME(/SECONDS )-t_loop
     
  ENDWHILE


  if keyword_set(print_times) then print,' Total time fminop : ',t_fminop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; POST-PROCESSING
  if keyword_set(print_times) then t=SYSTIME(/SECONDS)
  
  factor = 1D/total(x)
; print,  'factor', factor
  norm_x = x*factor
  achix = (reform(H#norm_x, operators.n_bases, n_elements(cmdata))) ; TF
  abs_hx = abs(achix)
  abs2_hx = abs2(achix, VERSION = version)
  arg_hx = angle(achix, VERSION = version)

  wx1 = abs2_hx*(cmdata.w_tan-cmdata.w_rad)
  wx2 = 2D*abs_hx*abs_y_data*cmdata.w_rad

  weights_x = {abs_hx:abs_hx, arg_hx:arg_hx, wx1:wx1, wx2:wx2}

  IF (display NE 0) THEN BEGIN
     wset, display+1
        WISARD_PLOT_FIT, RAD_FREQS = rad_freqs,$ 
                         ABS_HX = reform(abs_hx, n_elements(abs_hx)), $
                         ABS_Y  = reform(weights_constant.abs_y_data, n_elements(abs_hx)), $
                         CLOT_FREQS = clot_freqs,$
                         CLOT_FROM_DATA = reform(data.clot, n_elements(data.clot)),$
                         CLOT_FROM_CMDATA = reform(cloture_from_cmdata, n_elements(cloture_from_cmdata)),$
                         CLOT_FROM_CURR_X = reform(cloture_from_current_x, n_elements(cloture_from_cmdata)),$
                         INV_SIGMA_RAD = reform(sqrt(cmdata.w_rad),n_elements(abs_hx)), $
                         INDX_BADFLAG=indx_badflag,$ ; uncomment if bad flagged data is to be plotted
                         VERSION = version
  ENDIF

  norm_x = reform(norm_x, prior.np, prior.np)
  aux_output= {fov: fov, nbiter:nbiter, guess:nguess, np_min:np_min, resolution:fov/np_min, cmdata:cmdata, prior:prior, crit_array:crit_array, $
                                              x:norm_x, alpha:alpha, $
                                               operators:operators, weights_constant:weights_constant, $
                                               weights_alpha:weights_alpha, weights_x:weights_x, $
                                               rad_freqs:rad_freqs}
  
  if keyword_set(print_times) then print,' Time: Post-processing : ',SYSTIME(/SECONDS )-t

  return, norm_x
END
