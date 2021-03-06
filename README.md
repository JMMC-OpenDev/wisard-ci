WISARD-CI: the simplified interface to WISARD.
==============================================

[WISARD](http://www.mariotti.fr/doc/approved/JMMC-MAN-2500-0001.pdf) is a monochromatic image reconstruction program
developed for JMMC by S. Meimon and L. Mugnier, available under the
CeCILL-B license (see LICENSE and See Licence_CeCILL-B_V1-en.txt).

WISARD-CI is composed of WISARD plus a serie of enhancements and
additionnal procedures by G. DUVERT, making WISARD conformant to the
[Optical Imaging interface requirements](https://github.com/emmt/OI-Imaging-JRA/raw/master/doc/interface/OI-Interface.pdf),
to the [OI-FITS format](https://arxiv.org/abs/1510.04556), and to the [OImaging GUI](https://www.jmmc.fr/english/tools/data-analysis/oimaging/) available at
[](https://www.jmmc.fr/english/tools/data-analysis/oimaging/)

WISARD-CI installation: 
-----------------------

1. get the SVN repository files
2. put bin/wisard-ci in the PATH:
      > PATH=$PATH:/where/is/wisard-ci; export PATH 
3. or make a soft link from a $PATH directory to it: 
      > ln -s /where/is/wisard-ci ~/bin
   (provided ~/bin exist and is in the PATH of course)

- if you have IDL, you are rich. And you do not have to follow the instructions below. Go to "USAGE".

- if you do not have IDL, you are probably smart and have installed GDL instead.

  - install GDL, the free IDL clone, available as 'gnudatalanguage' in
  your ditribution (MacOS: "brew tap Homebrew/homebrew-science" and
  "brew install gnudatalanguage" for example; look for
  "gnudatalanguage" in debian, ubuntu etc.)

  - eventually, [get the source distribution of GDL](https://github.com/gnudatalanguage/gdl) and follow
  compilation instructions.

USAGE:
-------

if installation is ok, type "wisard-ci" and read help.

At that moment, you can either use wisard-ci in command line, or tell
the GUI that you have a local version of wisard-ci.

You can also call GDL (or IDL) and at the prompt, use the procedure wisardgui.pro.
use 
GDL> cd,"where/is/wisard-ci"; doc_library,"wisardgui" 
to print the minimum information on wisardgui procedure.

If you have IDL installed: It is the same as for GDL, but you must
edit bin/wisard-ci to replace "gdl" by "idl" and possibly add the
necessary IDL environment variables (not tested).

To USE a LOCAL INSTALL of WISARD-CI in the OImaging Gui:
--------------------------------------------------------

This is possible, provided you update the !PATH value in
wisard-ci/gdl_startup.pro : it is necessary to find the IDL procedures
of idlastro distribution, available at https://idlastro.gsfc.nasa.gov/

You start the local wisard-ci by selecting "WISARD" and not "WISARD (remote)" in the OImaging GUI.  


WARNING!: 

1) WISARD makes use of the astrolib (https://idlastro.gsfc.nasa.gov/)
set of procedures so be sure to have them in your (idl or gdl) !PATH.


2) WISARD uses precompiled OPTIMPACK libraries. See
https://cral.univ-lyon1.fr/labo/perso/eric.thiebaut/?Software/OptimPack
for details. If you have not a compatible library you'll need to build
one. And you are on your own on this.
