cd ..
rm -rf nexus-tools_build
mkdir nexus-tools_build
cd nexus-tools_build
export CXX=/usr/bin/g++-5
cmake -DLINUX_STATIC:BOOL=OFF -DCMAKE_BUILD_TYPE=Release ../nexus-tools
make flexcat
make nexcat
make MappingAnalyzer
make 5PrimeEndCounter
cp -f bin/flexcat ../nexus-tools/bin/
cp -f bin/nexcat ../nexus-tools/bin/
cp -f bin/MappingAnalyzer ../nexus-tools/bin/
cp -f bin/5PrimeEndCounter ../nexus-tools/bin/


