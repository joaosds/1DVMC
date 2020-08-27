#!/bin/bash

gfortran -O3 vmc1dHS.f90 -llapack -o vmc1dHS.exe
sleep 2 # wait 2 seconds
nohup time ./vmc1dHS.exe & 
