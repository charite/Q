cd ..
rm -rf preprocessing_build
mkdir preprocessing_build
cd preprocessing_build
export CXX=/usr/bin/g++-5
cmake -DLINUX_STATIC:BOOL=ON -DCMAKE_BUILD_TYPE=Release -DAVX2:BOOL=ON ../preprocessing
make flexcat
make nexcat
cp -f bin/flexcat ../../bin/
cp -f bin/nexcat ../../bin/



