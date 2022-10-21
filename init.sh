#!/bin/bash

DISPLAY=
#DISPLAY=-DDISPLAY

cd ./sequential/
make clean
make DISPLAY="$DISPLAY"

cd ../OpenMP/
make clean
make DISPLAY="$DISPLAY"

cd ../MPI/
make clean
make DISPLAY="$DISPLAY"

cd ../MPI_OpenMP/
make clean
make DISPLAY="$DISPLAY"
