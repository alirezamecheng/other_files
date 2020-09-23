#!/bin/sh
echo $JULIA_NUM_THREADS;
echo "Inter the julia thread number:";
int n;
read n;
export JULIA_NUM_THREADS=$n;
echo $JULIA_NUM_THREADS;
# julia1.1.0 example.jl;
