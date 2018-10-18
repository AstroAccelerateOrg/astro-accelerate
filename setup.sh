#!/usr/bin/env bash

echo "Setting up Astro-Accelerate..."

# CUDA required PATH and LD_LIBRARY_PATH
AA_ARCHITECTURE="64" #Leave blank for 32-bit, set to 64 for 64-bit
AA_ADD_PATH=/usr/local/cuda/bin
AA_ADD_LD_LIBRARY_PATH=/usr/local/cuda/lib${AA_ARCHITECTURE}

# Check if the proposed environment variable path exists
if [ ! -d "${AA_ADD_PATH}" ]; then
    echo "* ERROR: Directory ${AA_ADD_PATH} does not exist. Setup exiting..."
else
    # Check if PATH environment variable has been set on this system
    if [[ -z "${PATH}" ]]; then
	echo "* NOTICE: PATH variable is not yet set."
	echo "* NOTICE: PATH will be set to ${AA_ADD_PATH}"
	export PATH=${AA_ADD_PATH}
    else
	# PATH variable already exists"
	if [ -d "$AA_ADD_PATH" ] && [[ ":$PATH:" != *":$AA_ADD_PATH:"* ]]; then
	    echo "* NOTICE: Adding CUDA path to PATH environment variable."
	    export PATH="${PATH:+"$PATH:"}$AA_ADD_PATH"
	elif [[ $PATH = *"${AA_ADD_PATH}"* ]]; then
	    echo "* NOTICE: PATH already contains CUDA path."
	else
	    echo "* ERROR: Could not add ${AA_ADD_PATH} to PATH environment variable."
	    echo "* Either the requested path could not be found, or it has already been added to your PATH environment."
            echo "* Please check the CUDA version and path on your system."
	fi
    fi
fi

# Check if the proposed environment variable path exists
if [ ! -d "${AA_ADD_LD_LIBRARY_PATH}" ]; then
    echo "* ERROR: Directory ${AA_ADD_LD_LIBRARY_PATH} does not exist. Setup exiting..."
else
    # Check if LD_LIBRARY_PATH environment variable has been set on this system
    if [[ -z "${LD_LIBRARY_PATH}" ]]; then
        echo "* NOTICE: LD_LIBRARY_PATH variable is not yet set."
        echo "* NOTICE: LD_LIBRARY_PATH will be set to ${AA_ADD_LD_LIBRARY_PATH}"
        export LD_LIBRARY_PATH=${AA_ADD_LD_LIBRARY_PATH}
    else
       	# LD_LIBRARY_PATH variable already exists"
        if [ -d "$AA_ADD_LD_LIBRARY_PATH" ] && [[ ":$LD_LIBRARY_PATH:" != *":$AA_ADD_LD_LIBRARY_PATH:"* ]]; then
            echo "* NOTICE: Adding CUDA library path to LD_LIBRARY_PATH environment variable."
            export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+"$LD_LIBRARY_PATH:"}$AA_ADD_LD_LIBRARY_PATH"
        elif [[ $LD_LIBRARY_PATH = *"${AA_ADD_LD_LIBRARY_PATH}"* ]]; then
            echo "* NOTICE: LD_LIBRARY_PATH already contains CUDA path."
        else
            echo "* ERROR: Could not add ${AA_ADD_LD_LIBRARY_PATH} to LD_LIBRARY_PATH environment variable."
            echo "* Either the requested path could not be found, or it has already been added to your LD_LIBRARY_PATH environment."
            echo "* Please check the CUDA version and library path on your system."
        fi
    fi
fi

export CUDA_INSTALL_PATH=/usr/local/cuda/

echo "Done."
