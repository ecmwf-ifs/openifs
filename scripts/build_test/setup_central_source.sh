#!/usr/bin/env bash
#
# (C) Copyright 2011- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# This script will execute the ifstest create step within 
# OpenIFS, to generate source directory for a given cycle (branch). 
# Once created, the source directory will be removed and ifstest can be 
# run to set-up the links to the central source. 
#

set -eu
set -o pipefail

## Functions
log_echo ()
{
    if [[ -n ${setup_central_source_log-} ]]; then
        builtin echo "$@" | tee -a "${setup_central_source_log}"
    else 
        builtin echo "$@"
    fi
}

usage_info()
{  
    log_echo "----------------------------------------------------------------"
    log_echo "usage_info in setup_central_source called because" "$1"
    log_echo "----------------------------------------------------------------"
    log_echo "setup_central_source.sh sets up a directory in a user defined location" 
    log_echo "and then copies the standard openifs bundle dependencies, i.e. "
    log_echo "ecbuild, eccodes, multio, fckit, eckit, metkit, fcb5, atlas and fcm"
    log_echo ""
    log_echo "setup_central_source.sh usage :"
    log_echo "setup_central_source.sh -c <cycle>"
    log_echo ""
    log_echo "where :"
    log_echo "-c cycle. Default is 48r1"
    log_echo "This scripts requires source oifs-config.edit_me.sh, which "
    log_echo "is found in the openifs package directory"
    log_echo "----------------------------------------------------------------"
    exit 1 
}


check_user_input() {
    while true; do
        read -rp "$1 (y/n): " choice
        case "${choice,,}" in
        y | yes)
            return 0
            ;;
        n | no)
            return 1
            ;;
        *)
            echo "Invalid input. Please enter 'y' or 'n'."
            ;;
        esac
    done
}

check_oifs_config() {

    log_echo "[INFO]: Begin check oifs-config.edit_me.sh input"

    if [[ -n "${OIFS_HOME-}" && -d "${OIFS_HOME}" ]]; then 
        #
        # Use OIFS_HOME as the variable to check that oifs-config has been sourced
        #
        log_echo "[INFO]: $OIFS_HOME/oifs-config.edit_me.sh has been sourced"
        log_echo "[INFO]: Start setup of a central source for OpenIFS" 
        log_echo "[INFO]: OIFS_HOME directory path exists and it is $OIFS_HOME"
        #
        # Then create string with oifs_central_src + cycle
        #
        OIFS_CENTRAL_SRC_CYCLE=$OIFS_CENTRAL_SRC/$OIFS_CYCLE
        #
        # Check it exsits and create, if necessary, the directory to store the bundle source 
        #
        if [[ ! -d "$OIFS_CENTRAL_SRC_CYCLE" ]]; then
            log_echo "[WARNING]: $OIFS_CENTRAL_SRC_CYCLE does not exist."
            #
            # Request user input to check that the user is really sure.
            #   
            if check_user_input "[QUESTION]: Do you want to create $OIFS_CENTRAL_SRC_CYCLE ?"; then
                if ! mkdir -p "$OIFS_CENTRAL_SRC_CYCLE" 2>/dev/null ; then
                    log_echo "[ERROR]: Failed to create $OIFS_CENTRAL_SRC_CYCLE. Check permissions."
                    exit 1
                else
                    log_echo "[INFO]: $OIFS_CENTRAL_SRC_CYCLE created."
                fi
            else
                log_echo "[ERROR]: User chose not to create $OIFS_CENTRAL_SRC_CYCLE." 
                log_echo "         but this is required for this to continue - EXITING"
                exit 1
            fi
        else 
            log_echo "[INFO]: $OIFS_CENTRAL_SRC_CYCLE already exists, this will be overwritten"
        fi  
    else 
        log_echo "[ERROR]: OIFS_HOME does not exist or OIFS_HOME is not a directory."
        log_echo "         Please check oifs-config.edit_me.sh in the OpenIFS directory"
        log_echo "         and then source, e.g. source ~/openifs-48r1/oifs-config.edit_me.sh "
        log_echo "         EXITING"
        exit 1 
    fi

    # Also check that the openifs-test.sh script exists and is executable 
    if [[ -f "$OIFS_TEST/openifs-test.sh" && -x "$OIFS_TEST/openifs-test.sh" ]]; then
        log_echo "[INFO]: $OIFS_TEST/openifs-test.sh exists and it is executable" 
    else 
        log_echo "[ERROR]: $OIFS_TEST/openifs-test.sh is either missing or not executable"
        log_echo "         Please check and re-run - EXITING"
        exit 1
    fi

}

copy_bundle_dirs() {

    log_echo "[INFO]: find all directories $OIFS_HOME/source and copy to"
        log_echo "         $OIFS_CENTRAL_SRC/$OIFS_CYCLE"

    bundle_src_dir="$OIFS_HOME"/source

    for dir in "$bundle_src_dir"/*; do
        [[ -d "$dir" && ! -L "$dir" ]] || continue
        
        src_dir_name=$(basename "$dir")
        # NOTE: may want to check for existing dir and copy to backup
        # if [[ -d "$OIFS_CENTRAL_SRC/$OIFS_CYCLE/$src_dir_name" ]]; then 
        #     log_echo "[INFO]: move $OIFS_CENTRAL_SRC/$OIFS_CYCLE/$src_dir_name to .BAK"
        #     mv "$OIFS_CENTRAL_SRC/$OIFS_CYCLE/$src_dir_name" "$OIFS_CENTRAL_SRC/$OIFS_CYCLE/$src_dir_name.BAK"
        # fi
        log_echo "[INFO]: copy $dir to central source location: $OIFS_CENTRAL_SRC/$OIFS_CYCLE/"
        
        cp -rf "$dir" "$OIFS_CENTRAL_SRC/$OIFS_CYCLE/"
    done

    log_echo "[INFO]: Copy of source directories to $OIFS_CENTRAL_SRC/$OIFS_CYCLE/ complete"
}

check_central_source_dirs() {

    log_echo "[INFO]: Compare contents of $OIFS_CENTRAL_SRC/$OIFS_CYCLE/ to expected packages"
    log_echo "        defined in $OIFS_BUNDLE_SRC_DIR in main"

    remove_oifs_src=0

    for defined_dir in "${OIFS_BUNDLE_SRC_DIR[@]}"; do

        if [[ -d "$OIFS_CENTRAL_SRC/$OIFS_CYCLE/$defined_dir" ]]; then 

            log_echo "[INFO]: $defined_dir is in $OIFS_CENTRAL_SRC/$OIFS_CYCLE/"

        else 

            log_echo "[WARNING]: $defined_dir does not exist $OIFS_CENTRAL_SRC/$OIFS_CYCLE/"
            log_echo "           This could be a problem, do not remove OIFS_HOME/source"
            remove_oifs_src=1

        fi
    done

    for dir in "$OIFS_CENTRAL_SRC/$OIFS_CYCLE"/*; do
        src_dir_name=$(basename "$dir")

        # Skip directories ending with .BAK
        [[ "$src_dir_name" == *.BAK ]] && continue


        # Check if src_dir_name is in OIFS_BUNDLE_SRC_DIR
        if [[ " ${OIFS_BUNDLE_SRC_DIR[*]} " != *" $src_dir_name "* ]]; then
            log_echo "[WARNING]: central source directory, $src_dir_name does not exist in defined dirs"
            log_echo "           This may be a problem, do not remove OIFS_HOME/source - Check"
            remove_oifs_src=1
        fi
    done

    log_echo "[INFO]: Source directories copied to $OIFS_CENTRAL_SRC/$OIFS_CYCLE/ and checked"

    return "$remove_oifs_src"
    
}

main() {

    setup_central_source_log="setup_central_source.log"

    # This is the list of libraries that most OpenIFS users will 
    # not change, hence they can be downloaded and stored in a central
    # location.
    # NOTE: this list is valid for 48r1 and 49r1 but may not be valid for 
    #       later cycles. In the script, the list is derived by moving all 
    #.      directories and then compared to this list.  
    OIFS_BUNDLE_SRC_DIR=(
        "atlas"
        "ecbuild"
        "eccodes" 
        "eckit" 
        "fckit" 
        "fcm" 
        "fdb5" 
        "metkit" 
        "multio"
        )    

    while getopts :h OPTION ; do
        case "${OPTION}" in 
            h) 
                usage_info "[INFO]: user help requested "
            ;;
            \?)
                usage_info "[ERROR]: Invalid argument option -$OPTARG"
            ;;
        esac
    done

    # Check that oifs-config has been sourced, since this script depends on 
    # some of the environment variables
    check_oifs_config

    # Following checks, now do the business...
    
    # cd into the openifs package directory
    cd "$OIFS_HOME"
    
    # run openifs-test create step to set up the sources for an openifs build
    log_echo "[INFO] : Running $OIFS_TEST/openifs-test.sh create step, which "
    log_echo "         will create $OIFS_HOME/source"
    
    "$OIFS_TEST"/openifs-test.sh -c

    # Identify all the directories to copy and copy them and then check against 
    # OIFS_BUNDLE_SRC_DIR
    copy_bundle_dirs 

    if check_central_source_dirs; then 
        log_echo "[INFO]: Directories copied and checked with no warnings, now remove $OIFS_HOME/source"
        if [[ -d ${OIFS_HOME}/source ]]; then 
            if rm -rf "${OIFS_HOME}/source"; then
                log_echo "[INFO]: Successfully removed directory: ${OIFS_HOME}/source"
            else
                log_echo "[ERROR]: Failed to remove directory: ${OIFS_HOME}/source"
            
            fi
        fi
    else
        log_echo "[INFO]: Directories copied but checks produced warnings, do not remove $OIFS_HOME/source"
        log_echo "        Please check $OIFS_HOME/source and $OIFS_CENTRAL_SRC/$OIFS_CYCLE/"
    fi
 
 }

main "$@"
