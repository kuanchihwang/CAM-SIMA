#! /bin/bash

# Script to run doctest tests and unit tests on CAM python scripts
# Note: This script should be kept up to date with
#       .github/workflows/pr_open_sync_workflow.yml

# We only test python3 now
PYTHON=python3
NUMTESTS=0
NUMERRORS=0

run_doctest() {
  # Run the doctests in the input python module ($1)
  echo "Running doctests in ${1}"
  if [ `grep -c 'import doctest' ${1}` -ne 0 ]; then
    ${PYTHON} ${1}
    res=$?
  else
    ${PYTHON} -m doctest ${1}
    res=$?
  fi
  if [ $res -ne 0 ]; then
    echo "ERROR: doctest for ${1} returned ${res}"
    NUMERRORS=$(( NUMERRORS + 1 ))
  fi
  NUMTESTS=$(( NUMTESTS + 1 ))
}

run_unittest() {
  # Run the unit tests in the input python module ($1)
  echo "Running unit tests in ${1}"
  ${PYTHON} ${1}
  res=$?
  if [ $res -ne 0 ]; then
    echo "ERROR: unittest for ${1} returned ${res}"
    NUMERRORS=$(( NUMERRORS + 1 ))
  fi
  NUMTESTS=$(( NUMTESTS + 1 ))
}

# Try to find the correct directory
scriptdir="$( cd $( dirname ${0} ); pwd -P )"
if [ ! -d "cime_config" ]; then
  # Script was not called from top level directory, try to reset
  cd $( dirname ${scriptdir} )
fi

# Make sure we are at the top level
if [ ! -d "cime_config" ]; then
  echo "ERROR: Cannot find CAM sandbox to test"
  exit 1
fi

# CAM config doctests:
run_doctest cime_config/cam_config.py
# CAM autogen doctests:
run_doctest cime_config/cam_autogen.py
# CAM build cache doctests:
run_doctest cime_config/cam_build_cache.py
# Registry generator doctests:
run_doctest src/data/generate_registry_data.py
# CAM config unit tests:
run_unittest test/unit/cam_config_unit_tests.py
# Registry generator unit tests:
run_unittest test/unit/test_registry.py
# Physics variable init (phys_init) generator unit tests:
run_unittest test/unit/write_init_unit_tests.py

# Report
if [ ${NUMERRORS} -gt 0 ]; then
  echo "${NUMERRORS} out of ${NUMTESTS} tests FAILED"
else
  echo "All ${NUMTESTS} tests PASSED!"
fi