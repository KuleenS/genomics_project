#!/bin/bash

cmake --build build

./build/alignment test_input.txt global 0 penalty.csv
./build/alignment test_input.txt global 1 penalty.csv
./build/alignment test_input.txt global 2 penalty.csv

./build/hirschberg test_input.txt penalty.csv
./build/four_russians test_input.txt 2