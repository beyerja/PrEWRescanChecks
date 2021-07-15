#!/bin/bash

echo "#########################################################################"
echo "### Loading needed software versions ####################################" 
echo "#########################################################################"

lcg_version=LCG_100
if [[ "${LCG_VERSION}" == "${lcg_version}" ]]; then 
  echo "LCG view already loaded."
elif [ ! -z "${LCG_VERSION}" ]; then
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "ERROR: Another LCG view is already loaded!"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
else 
  echo "Loading LCG view - version: ${lcg_version}" 
  source /cvmfs/sft.cern.ch/lcg/views/${lcg_version}/x86_64-centos7-gcc10-opt/setup.sh
fi  