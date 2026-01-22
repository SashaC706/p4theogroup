cls
@REM g++ ./main.cpp -O3 -ffast-math -march=native -pthread --std=c++26 -o ./o.exe
@REM g++ ./RK4.cpp -O3 --std=c++26 -o ./o.exe
@REM g++ ./RK4.cpp --std=c++26 -o ./o.exe
g++ ./main.cpp ./source/constants.hpp ./source/LE.cpp ./source/project.hpp -o ./o.exe
o.exe