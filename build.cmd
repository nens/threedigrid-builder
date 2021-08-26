cmake -G "MSYS Makefiles" -DCMAKE_GNUtoMS=ON -DCMAKE_BUILD_TYPE="RELEASE" -B.\build .
cmake --build .\build --target install/strip
