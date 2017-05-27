#!/bin/bash
#echo "creating a xcode project from dcshell codes"
echo "== start =="

#cmake build_type_option code_path
mkdir -p build/debug
cd build/debug
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-g" ../..
cd -

mkdir -p build/release
cd build/release
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3" ../..
cd -

if [[ "$OSTYPE" == "darwin"* ]]; then
	mkdir -p build/xcode
	cd build/xcode
	# -G <generator name> <path to source>
	cmake -G Xcode ../..
	cd -
fi

cd build/release
#make -j4
make
cd -

ln -s build/release/dcshell .
echo "=== done ==="
