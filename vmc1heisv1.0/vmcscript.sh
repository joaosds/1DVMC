#!/bin/bash

gfortran -O3 vmc1dheis.f90 -llapack -o vmc1dheis.exe
sleep 2 # wait 2 seconds
nohup time ./vmc1dheis.exe & 
