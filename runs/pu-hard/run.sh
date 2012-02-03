#!/bin/bash


rm -f *.ini in.* coarse* log-* pair.* 
../../cg-ibi.py --block "HH" --num-blocks 1 --num-chains 500
