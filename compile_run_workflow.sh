#!/bin/bash

./gen_codes.sh

cd src; make

cd ..

cd runs

run.sh


