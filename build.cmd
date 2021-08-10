cmake -G "MSYS Makefiles" -A Win32 -DCMAKE_GNUtoMS=ON -DCMAKE_BUILD_TYPE="RELEASE" -B.\build .
cmake --build .\build --target install
