#!/bin/bash

# ==============================================================================
# Define usage
#     Install conda (if needed) and create the conda environment
# ==============================================================================
## define usage
function usage {
    echo "----------------------------------------------------------------"
    echo "Runs install scripts in install.d and sets the "
    echo "environment with an environment.d file"
    echo "----------------------------------------------------------------"
    echo ""
    echo "  FINITEFAULT_DIR string  the directory where the finite fault"
    echo "                          repository exists"
    echo "                          (example: /home/user/neic-finitefault)"
}

## parse arguments
REQUIRED_ARGS=()
while [[ $# -gt 0 ]]; do
    # shellcheck disable=SC2221,SC2222
    case $1 in
        -h|--help)
            usage
            exit 0
            ;;
        -*|--*)
            echo "Invalid option $1"
            usage
            exit 1
            ;;
        *)
            REQUIRED_ARGS+=("$1")
            shift
            ;;
    esac
done

### set positional arguments
set -- "${REQUIRED_ARGS[@]}" # restore positional parameters
FINITEFAULT_DIR=$1
if [ -z ${1+x} ]
then 
    echo "Argument FINITEFAULT_DIR must be set";
    exit 1;
fi


# Create finite fault env
if conda env list | grep "ff-env"; then
    echo "The miniforge environment 'ff-env' already exists. Delete and recreate it? [y/n]"
    # shellcheck disable=SC2162
    read response;
    if [ "$response" = "y" ] || [ "$response" = "Y" ]; then
        conda env remove -n ff-env -y;
        # install ff-env
        echo "Creating a finite fault environment called 'ff-env'"
        conda env create -f "${FINITEFAULT_DIR}"/install.d/environment.yml
        echo "Environment created. Use 'conda activate ff-env' to start the environment"
    elif [ "$response" = "n" ] || [ "$response" = "N" ]; then
        echo "Skipping environment install";
    else
        echo "Unknown response. Exiting.";
        exit 1;
    fi
else
    # install ff-env
    echo "Creating a finite fault environment called 'ff-env'"
    conda env create -f "${FINITEFAULT_DIR}"/install.d/environment.yml
    echo "Environment created. Use 'conda activate ff-env' to start the environment"
fi
