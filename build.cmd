cmake -G "MSYS2 Makefiles" -DCMAKE_GNUtoMS=ON -DCMAKE_BUILD_TYPE="RELEASE" -B.\build .
cmake --build .\build --target install
