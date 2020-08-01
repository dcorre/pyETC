#!/bin/bash
#
# Install pyETC and dependencies using Conda

echo "$(uname)";
# install miniconda
if ! which conda 1> /dev/null; then
    if test ! -f ${HOME}/miniconda/etc/profile.d/conda.sh; then
        # install conda
        [ "$(uname)" == "Darwin" ] && MC_OSNAME="MacOSX" || MC_OSNAME="Linux"
        MINICONDA="Miniconda${PYTHON_VERSION%%.*}-latest-${MC_OSNAME}-x86_64.sh"
        curl -L https://repo.continuum.io/miniconda/${MINICONDA} -o miniconda.sh
        #bash miniconda.sh -b -u -p ${HOME}/miniconda
        bash miniconda.sh -b -p ${HOME}/miniconda
    fi
    source ${HOME}/miniconda/etc/profile.d/conda.sh
fi
hash -r

# get CONDA base path
CONDA_PATH=$(conda info --base)

# configure miniconda
conda config --set always_yes yes --set changeps1 no
conda config --add channels conda-forge
conda update --quiet --yes conda
conda info --all

# create environment for tests (if needed)
if [ ! -f ${CONDA_PATH}/envs/pyETC/conda-meta/history ]; then
    conda create --name pyETC python=${PYTHON_VERSION} pip setuptools
fi

# install conda dependencies (based on pip requirements file)
conda run --name pyETC \
conda install --name pyETC --quiet --yes --file requirements_dev.txt --update-all

# install other conda packages that aren't represented in the requirements file
#conda install --name gwpyci --quiet --yes \
#    "python-framel>=8.40.1" \
#    "python-lal" \
#    "python-lalframe" \
#    "python-lalsimulation" \
#    "python-ldas-tools-framecpp" \
#    "python-nds2-client" \
#    "root>=6.20" \
#    "root_numpy" \
#;

# activate the environment
if MC_OSNAME="MacOSX"; then
	. ${CONDA_PATH}/etc/profile.d/conda.sh
fi
conda activate pyETC
