# EventsGenerator

This generator implements ideas developed at LAL by Emi Kou, Andrey Tayduganov and internship made by Andrii Kotenko in fall of 2017. It is about generation of 4-body radiative B-meson decay.

Prerequirities:

ROOT v5
gcc 4.8+
CMake 2.6+

How to run:

cd Generator
mkdir build
cd build
cmake ../
make -jN , where N - number of threads available
./main
