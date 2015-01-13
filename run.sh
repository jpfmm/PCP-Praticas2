#!/bin/sh

echo "Compiling..."

make

echo "Running the tests"

./bin/lab2

echo "You can now copy timings.dat"