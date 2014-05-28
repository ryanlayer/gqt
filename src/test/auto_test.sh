#!/bin/bash

while [ 1 ]
do
    if [ "test_genotq" -ot "test_genotq.c" ]
    then
        make clean
        make
    fi
    sleep 0.4
done
