#!/bin/bash

./raw 0.5 result500/ 100000000
./detected 0.5 result500/ 100000000
./smeared 0.5 result500/ 100000000
./cut 0.5 result500/ 100000000

./raw 0.25 result250/ 100000000
./detected 0.25 result250/ 100000000
./smeared 0.25 result250/ 100000000
./cut 0.25 result250/ 100000000

./raw 0.1 result100/ 100000000
./detected 0.1 result100/ 100000000
./smeared 0.1 result100/ 100000000
./cut 0.1 result100/ 100000000

./raw 0.05 result50/ 100000000
./detected 0.05 result50/ 100000000
./smeared 0.05 result50/ 100000000
./cut 0.05 result50/ 100000000

