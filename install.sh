#!/bin/sh
chmod +x ./*.py
cd symmetry_source
qmake Symmetry.pro -r -spec linux-g++
make -r -w in ../Symmetry_Release
cp Symmetry_Release/Symmetry ./
chmod +x Symmetry