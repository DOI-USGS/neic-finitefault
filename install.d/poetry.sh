#!/bin/bash

# ==============================================================================
# Define usage
#     
# ==============================================================================
## define usage
function usage {
    echo "----------------------------------------------------------------"
    echo "Install pip installable dependencies with poetry. "
    echo "NOTE: ff-env must be created with conda.sh first."
    echo "----------------------------------------------------------------"
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

### activate ff-env
conda activate ff-env;
conda env list

### install the main dependencies
poetry install;
source $(poetry env info --path)/bin/activate;

### install okada separately since it depends on numpy
poetry run poe okada;