; $Id: op_init.pro,v 1.5 2011/11/24 10:28:14 lmugnier Exp $
PRO OP_INIT, name, $
             prefix_output = prefix, $
             basename_output = basename, $
             suffix_output = suffix
;+
; NAME:
;   OP_INIT
;
;
; PURPOSE:
;   Initialize OptimPack routines.
;
;
; CATEGORY:
;   OptimPack.
;
;
; CALLING SEQUENCE:
;   OP_INIT [, name] [, PREFIX_OUTPUT = prefix]
;                    [, BASENAME_OUTPUT = basename]
;                    [, SUFFIX_OUTPUT = suffix]
;
;
; INPUTS:
;   NAME - Name (without extension) of the OptimPack-IDL library. If not
;          given, it defaults to "OptimPack_IDL${OSTYPE}" or
;          "OptimPack_IDL${OSTYPE}_64" for 64 bit platforms (currently Linux
;          MacOSX and Solaris), where ${OSTYPE} is computed by op_init. Its
;          value is "solaris" for Solaris, "linux" for Linux, "darwin" for Mac
;          OS X and "windows" for Windows.
;
;
; OPTIONAL INPUTS:
;   None.
;
; KEYWORD PARAMETERS:
;   PREFIX_OUTPUT   : (optional output) variable name that receives the prefix of the 
;                     library name on output ('' or './') 
;   BASENAME_OUTPUT : (optional output) variable name that receives the basename of the
;                     library name on output (name or or 'OptimPack_IDL'+OSTYPE) 
;   SUFFIX_OUTPUT   : (optional output) variable name that receives the suffix of the
;                     library name on output ('.so' or '.DLL' or etc...)
;
; OUTPUTS:
;   None.
;
; OPTIONAL OUTPUTS:
;   None.
;
;
; COMMON BLOCKS:
;   OP_COMMON - used to store the name of the shared OptimPack library.
;
;
; SIDE EFFECTS:
;   Data in common block OP_COMMON is updated.
;
;
; RESTRICTIONS:
;   None.
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;   op_init, "OptimPack_IDL"
;
;
; MODIFICATION HISTORY:
;   $Log: op_init.pro,v $
;   Revision 1.5  2011/11/24 10:28:14  lmugnier
;   Added compatibility with 32bit-IDL on Mac OS X.
;   For Mac OS X library name is thus either "OptimPack_IDLdarwin.dylib" or
;   "OptimPack_IDLdarwin_64.dylib".
;
;   Revision 1.4  2010/05/05 16:40:06  lmugnier
;   Fixed Mac OS X compatibility (thanks to Xavier Haubois) and updated
;   documentation. Library name is "OptimPack_IDLdarwin_64.dylib" for Mac OS X.
;
;   Revision 1.3  2008/10/24 16:50:20  mugnier
;   - Added PREFIX_OUTPUT, BASENAME_OUTPUT, SUFFIX_OUTPUT to output these variables.
;   - Added windows and MacOSX compatibility.
;   - Removed VMS, no longer supported by IDL.
;
;   Revision 1.2  2006/04/26 15:36:02  mugnier
;   - Major change: changed default library basename from 'OptimPack_IDL' to
;    'OptimPack_IDL${OSTYPE}' for multi-platform use.
;
;   - Added support for Linux 64 bits platforms.
;
;   Revision 1.1  2006/04/26 15:28:20  mugnier
;   Initial revision
;
;   Initial version : 2003, Eric THIEBAUT.
;   put under version control as 1.1 (see Log).
;-
  common op_common, libname
  on_error, 2
  
;  Notes: 
; getenv('OSTYPE')='' under windows so don't use getenv
; OSTYPE = STRLOWCASE(!version.OS_NAME) yields 2 words ("microsoft windows") so don't use it
; OSTYPE = STRLOWCASE(!version.OS) yields "sunos" for "solaris" so don't use it like this
  prefix = ''
  CASE !version.os_family OF
     'unix': BEGIN
        OSTYPE = STRLOWCASE(!version.OS_NAME) ; solaris for sunos, linux for linux,
                              ; 'Mac OS X' for MacOSX.
        OS = STRLOWCASE(!version.os)

      CASE OS OF
        'hp-ux': suffix = '.sl'
        'aix': BEGIN
          ;; AIX won't find a shared lib in the current dir
          ;; unless the name is preceded with a ./
          suffix = '.a'
          IF strmid(basename, 0, 1) ne '/' THEN prefix = './'
       END 
        'sunos': BEGIN
          IF (!version.memory_bits eq 64) THEN suffix = '_64.so' $
          ELSE                                 suffix = '.so'
       END 
        'linux': BEGIN
          IF (!version.memory_bits eq 64) THEN suffix = '_64.so' $
          ELSE                                 suffix = '.so'
       END
        'darwin': BEGIN
          IF (!version.memory_bits eq 64) THEN suffix = '_64.dylib' $ ; for DYnamic LIBrary.
          ELSE                                 suffix = '.dylib'
          OSTYPE = 'darwin'  ; 
        END
        ELSE: suffix = '.so'
     ENDCASE
   END 
;    'vms':     suffix = '.EXE' ; VMS support dropped w/ IDL 5.5
    'Windows': BEGIN 
       suffix = '.DLL'
       OSTYPE = 'windows'
    END
    ELSE: message, "Don't know what to do with: " + !version.os_family
  ENDCASE
  
  IF n_elements(name) eq 0L THEN BEGIN
    basename = 'OptimPack_IDL'+OSTYPE
  ENDIF  ELSE BEGIN
    basename = name
  ENDELSE 
  libname = prefix + basename + suffix
END 
