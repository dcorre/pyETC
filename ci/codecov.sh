#!/bin/bash

set -ex
trap 'set +ex' RETURN

#
# Submit coverage data to codecov.io
#

# reactivate environmennt
if [ -n ${CIRCLECI} ] && [ -d /opt/conda/envs ]; then
    conda activate test-env || source activate test-env 
fi

# get path to python
PYTHON_VERSION=$(echo "${PYTHON_VERSION:-${TRAVIS_PYTHON_VERSION}}" | cut -d. -f-2)
PYTHON=$(which "python${PYTHON_VERSION}")

# install codecov
${PYTHON} -m pip install ${PIP_FLAGS} coverage codecov

# find job name
_JOBNAME=${CIRCLE_JOB:-${TRAVIS_JOB_NAME}}

# submit coverage results
${PYTHON} -m codecov --file tests/coverage.xml --flags $(uname) python${PYTHON_VERSION/./} ${_JOBNAME%%:*}
