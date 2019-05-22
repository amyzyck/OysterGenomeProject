#! /bin/bash

for d in */ ; do
    cd "$d"
    convert manhattan_*.png manhattan_plots_combined.pdf
    convert outlier_*.png outlier_compare_combined.pdf
    cd ..
done