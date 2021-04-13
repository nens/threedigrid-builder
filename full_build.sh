#!/bin/bash

set -e
if [ $# -eq 0 ]
then  
    BUILD_TYPE="RELEASE"
else
	BUILD_TYPE=${1^^}
	echo $BUILD_TYPE
fi

cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE . -B/tmp/build
cmake --build /tmp/build --clean-first --target install
