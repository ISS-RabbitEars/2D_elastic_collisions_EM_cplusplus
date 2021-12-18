#!/bin/bash

mkdir pos
g++ collisions.cpp -o col
./col
rm col
cp render.sh pos/
cp pos_template.pov pos/
cd pos
./render.sh
#ffmpeg -framerate 60 -i frame_%d.png 2dcol_pbc.gif
ffmpeg -framerate 60 -i frame_%d.png colsim.mp4
rm render.sh
rm pos_template.pov
rm N.dat
rm rs.dat
rm q.dat
rm nframes.dat
