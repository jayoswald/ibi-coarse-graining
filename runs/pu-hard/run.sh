#!/bin/bash


rm -f *.ini in.* coarse* log-* pair.*  *.lammpstrj restart*
../../cg-ibi.py --block "HH" --num-blocks 1 --num-chains 500 --system_type hard
