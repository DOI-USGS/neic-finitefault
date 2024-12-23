#!/bin/bash

# ==============================================================================
# Define usage
#     Install miniforge (if needed)
# ==============================================================================
## define usage
function usage {
    echo "----------------------------------------------------------------"
    echo "Install Conda"
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


# shellcheck source=./install.d/miniforge.sh
source ./install.d/miniforge.sh;
