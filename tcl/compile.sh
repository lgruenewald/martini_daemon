#!/bin/sh

gcc -shared $(pkg-config --cflags tcl8.6) -o bond_loader.so bond_loader.c $(pkg-config --libs tcl8.6) -fPIC -g
