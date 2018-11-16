#!/bin/bash
gcc -g -c math_functions.c -o math_functions.o
gcc -g -c separation_v0.c -o separation.o
gcc -g SOBI.c -lm math_functions.o separation.o -o ./SOBI.exe
./SOBI.exe
