cmake -G"MinGW Makefiles" -DCMAKE_BUILD_TYPE="RELEASE" -B.\build .
cmake --build .\build --target install
