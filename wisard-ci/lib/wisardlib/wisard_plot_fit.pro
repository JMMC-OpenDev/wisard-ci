; $Id: wisard_plot_fit.pro,v 1.11 2010-12-08 14:56:08 duvert Exp $ 
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
;   edited by J�r�me Idier, ISTE, London, 2008.
;

;+
; NAME:
;	WISARD_PLOT_FIT - plots the fitting between data and reconstruction in visibility moduli 
;
; CATEGORY:
;	Array and Image Processing Routines 
;
; CALLING SEQUENCE:
;	WISARD_PLOT_FIT,RAD_FREQS = rad_freqs, ABS_HX = abs_hx, ABS_Y = abs_y, $
;             INV_SIGMA_RAD = inv_sigma_rad,INDX_BADFLAG=indx_badflag,
;             COLORS = colors [, /VERSION] [, /HELP]
;
;
; PURPOSE:
;   This procedure plots the visibility moduli for the data (abs_y), the model (abs_hx), and the
;residual in moduli (|abs_y-abs_hx|) in tenths of sigma.
;
;
; POSITIONAL PARAMETERS:
;   None.
;
; KEYWORD PARAMETERS:
;
;	RAD_FREQS     : (input) radial frequencies
;
;	ABS_HX        : (input) visibility moduli for the reconstruction
;
;	ABS_Y         : (input) visibility moduli for the (myopic convexified) data
;
;	INV_SIGMA_RAD : (input) inverse of the radial standard deviation for the
;                   myopic convexified data model
;                   
;  INDX_BADFLAG : (optional input) vector containing the indices of "bad flag" data. 
;                 If present, bad flag data is plotted in a different color, otherwise
;                 it is not plotted (default).
;
;   COLORS   : (optional input) 3-element vector containing the "grey levels"
;              of the plots. Default = [254, 150, 100], appropriate for use
;              with rainbow+white color table (39). The 1st element is used for
;              reconstructed visibilities, the 2nd element for measured
;              visibilities, the 3rd element for the difference.
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
;   edited by J�r�me Idier, ISTE, London, 2008.
;
; EXAMPLE:
;when wisard has finished and if aux_output exists, you can draw a modulus
;visibility fitting plot with:
;
;rad_freqs=aux_output.rad_freqs
;abs_hx=aux_output.weights_x.abs_hx
;abs_hx=reform(abs_hx,n_elements(abs_hx))
;abs_y=aux_output.weights_constant.abs_y_data
;abs_y=reform(abs_y,n_elements(abs_y))
;inv_sigma_rad = reform(sqrt(aux_output.cmdata.w_rad),n_elements(abs_hx))
;WISARD_PLOT_FIT,RAD_FREQS = rad_freqs, ABS_HX = abs_hx, ABS_Y = abs_y, $
;             INV_SIGMA_RAD = inv_sigma_rad, VERSION = version, HELP = help
             
; SEE ALSO:
;   All the WISARD_*.pro files, which are part of WISARD and covered by the
;   same copyright and license.
;   
; HISTORY:
;   $Log: not supported by cvs2svn $
;   Revision 1.10  2010/10/25 16:56:50  mvannier
;   Commit general au niveau wisard/
;
;  
;   Revision 1.91 2010/03/16 mvannier
;   INDX_BAD_FLAG is now an optional input. If not set, bad-flagged data are not plotted
;   
;   Revision 1.9  2009/11/24 10:50  vannier
;   INDX_BAD_FLAG added as input and used for different display of bad-flagged observed and 
;   reconstructed data.
;
;   Revision 1.8  2008/01/07 16:56:26  mugnier
;   Added COLORS optional input keyword.
;
;   Revision 1.7  2007/10/31 12:44:20  mugnier
;   Added 3rd reference.
;
;   Revision 1.6  2007/02/01 13:31:31  meimon
;   *** empty log message ***
;
;   Revision 1.5  2007/01/26 19:08:10  meimon
;   doc added
;
;-

PRO WISARD_PLOT_FIT,RAD_FREQS = rad_freqs, ABS_HX = abs_hx, ABS_Y = abs_y, $
                    CLOT_FREQS = CLOT_FREQS, $
                    CLOT_FROM_DATA =CLOT_FROM_DATA, $
                    CLOT_FROM_CMDATA =CLOT_FROM_CMDATA, $
                    CLOT_FROM_CURR_X =CLOT_FROM_CURR_X, $
                    INV_SIGMA_RAD = inv_sigma_rad, INDX_BADFLAG=indx_badflag, $
                    COLORS = colors, $
                    VERSION = version, HELP = help

;  on_error,2
  T3DATA=CLOT_FROM_DATA*180./!DPI
  T3CMDATA=CLOT_FROM_CMDATA*180./!DPI
  T3XDATA=CLOT_FROM_CURR_X*180./!DPI

  IF keyword_set(version) THEN $
     printf, -2, '% '+routine_courante()+': $Revision: 1.11 $, $Date: 2010-12-08 14:56:08 $'

  IF (n_params() NE 0) OR keyword_set(help) THEN BEGIN ;(fill-in from syntax)
     message, 'Help required or incorrect syntax. Documentation:', /INFO
     doc_library, routine_courante()
     retall
  ENDIF

  xtitle =  'frequency'         ; textoidl('|\nu|') 
  ytitle = 'VIS.'
  xmax=max(rad_freqs)*1.05D
  xmin=0.0D
  IF (n_elements(colors) EQ 0L) THEN $
     colors_inside = [254, 150, 100, 60] $
  ELSE $
     colors_inside = colors
  !P.MULTI=[0,0,2,0,0]
  !P.POSITION=[0.05,0.3,0.95,0.95]
; flags. my position is to plot only good data.
  indx_goodflag=indgen(n_elements(rad_freqs))
  if (keyword_set(indx_badflag)) then if (min(indx_badflag) GE 0) then begin
     indx_goodflag[indx_badflag]=-1
     w=where(indx_goodflag eq -1, comp=indx_goodflag)
  end

  plot,rad_freqs[indx_goodflag],abs_y[indx_goodflag],  XRANGE = [xmin, xmax], XMARGIN=[10,4], $
       XTICKS=1,  XSTYLE=1, color = -1  ,/NODATA ;, /YLOG ; white, box only

  oplot,rad_freqs[indx_goodflag],abs_y[indx_goodflag],psym=6,color = colors_inside[1] ;green
  oplot,rad_freqs[indx_goodflag],abs_hx[indx_goodflag],psym=7, color = colors_inside[0] ; red

  tstring = ['Abs(Measured Vis.)','Abs(Reconstructed Vis.)']
  tlinestyle= [0,0]
  tsym= [7, 6]
  if (keyword_set(indx_badflag)) then if (min(indx_badflag) GE 0) then begin
     tstring = ['Abs(Reconstructed Vis.)', 'Abs(Measured Vis.)','Flagged data']
     tlinestyle= [0,0,0]
     tsym= [7, 6, 6]
     oplot,rad_freqs[indx_badflag],abs_y[indx_badflag],psym=6,color = colors_inside[3] ;blue
     oplot, rad_freqs[indx_badflag],abs_hx[indx_badflag],psym=7, color = colors_inside[3] ;blue
  endif

  !P.COLOR=-1
  wis_legend, tstring, line = tlinestyle, psym = tsym, color = colors_inside, /bottom, /left, clear = clear

  !P.POSITION=[0.05,0.1,0.95,0.3]
  plot,clot_freqs,T3DATA, color=-1,  XRANGE = [0, max(clot_freqs)*1.05], /NODATA
  oplot,clot_freqs,T3DATA, color=colors_inside[1] , psym=6 ; green
  oplot, clot_freqs, T3XDATA,psym=7,color =colors_inside[0] ; red
;  oplot, clot_freqs, T3CMDATA, color=colors_inside[1], psym=5
  tstring = ['Measured closures','Reconstructed closures'];, 'Closure from initial cmdata',]
  tlinestyle= [0,0];,0]           ;
  tsym= [7, 6];, 5]
  wis_legend, tstring, line = tlinestyle, psym = tsym, color = colors_inside, /bottom, /left, clear = clear

;  ord_rad = sort(rad_freqs)
;  plot, rad_freqs,(1.0D*abs(abs_hx-abs_y)*inv_sigma_rad)[ord_rad], XRANGE = [xmin, xmax], XMARGIN=[10,4], $
;        XTIT = xtitle, XSTYLE=1, color = -1 ,/NODATA
;  oplot, rad_freqs[ord_rad],(1.0D*abs(abs_hx-abs_y)*inv_sigma_rad)[ord_rad], color=colors_inside[2] ;light blue
;  tstring = ['Abs(Difference)/standard deviation']
;  !P.COLOR=-1
;  tlinestyle= [0]
;  tsym= [-3]
;  wis_legend, tstring, line = tlinestyle, psym = tsym, color = colors_inside[2], /top, $
;              /right, clear = clear

  !P.MULTI=0
  !P.POSiTION=0
END
