#!/bin/bash
rm -f main

for i in {1..10}; do echo ' '; done

make all

for i in {1..10}; do echo ' '; done

./main
