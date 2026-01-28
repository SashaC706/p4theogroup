clear
# ++ ./main.cpp -O3 -ffast-math -march=native -pthread --std=c++26 -o ./o.exe
g++ ./main.cpp ./source/constants.hpp ./source/LE.cpp ./source/project.hpp -march=native -pthread --std=c++26 -o ./o.exe
./o.exe