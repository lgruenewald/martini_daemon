#!/bin/sh

gcc -shared $(pkg-config --cflags tcl) -o bond_loader.so bond_loader.c $(pkg-config --libs tcl) -fPIC -g
