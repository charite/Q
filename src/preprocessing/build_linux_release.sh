cd ..
rm -rf preprocessing_build
mkdir preprocessing_build
cd preprocessing_build
export CXX=/usr/bin/g++-5
cmake -DLINUX_STATIC:BOOL=OFF -DCMAKE_BUILD_TYPE=Release ../preprocessing
make flexcat
make nexcat
mkdir ../../bin
cp -f bin/flexcat ../../bin/
cp -f bin/nexcat ../../bin/


