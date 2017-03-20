#!/bin/bash
function printUsage ()
{
    echo -e "Usage: $0 [-h] [-d] [-f fov] [-i init_img] oifits_input oifits_output"
    echo -e "          -h : help"
    echo -e "          -d : display graphs"
    echo -e "          -i init_img : pass startup (guess) image (FITS format)"
    echo -e "          -f fov : field-of-view in mas"
    echo -e "          -n nbiter : number of iterations"
    exit 1
}

# Parse command-line parameters
while getopts "hdf:i:n:" option
do
    case $option in
        h ) # Help option
            printUsage;
            exit 0;;
        f ) # fov
            FOVCOMMAND=', fov='"$OPTARG"' ';
	    ;;
        d ) # display
            DISPLAYCOMMAND=', /display ';
	    ;;
        i ) # init_image
            INIT_IMAGE_COMMAND=', init_img='\'"$OPTARG"\'' ';
	    ;;
        n ) # nbiter
            NBITER_COMMAND=', nbiter='\'"$OPTARG"\'' ';
	    ;;
        * ) # Unknown option
            echo "Invalid option -- $option"
            printUsage;
            exit 1;;
    esac
done
let SHIFTOPTIND=$OPTIND-1
shift $SHIFTOPTIND

export LD_LIBRARY_PATH=/opt/idl/idl64/bin/bin.linux.x86_64::/home/gildas/INTROOT/lib
export IDL_DIR=/opt/idl/idl64
export IDL_STARTUP=/home/gildas/amber/idl/idl_startup.pro

cd $HOME/wisard-ci
LIGNE="wisardgui,'"$1"','"$2"'"${DISPLAYCOMMAND}${FOVCOMMAND}${INIT_IMAGE_COMMAND}${NBITER_COMMAND}
# /opt/idl/idl64/bin/idl -e $LIGNE
#echo $LIGNE
 gdl -e "$LIGNE"

