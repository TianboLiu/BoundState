#!/bin/bash
Eq="$1"

./virtual ${Eq} 0
./virtual ${Eq} 0.001
./virtual ${Eq} 0.005
./virtual ${Eq} 0.01
./virtual ${Eq} 0.05
./virtual ${Eq} 0.1
./virtual ${Eq} 0.5

unset Eq
