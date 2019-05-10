#! /bin/bash

for d in */ ; do
    cd "$d"
    convert *.png manhattan_plots_combined.pdf
    cd ..
done