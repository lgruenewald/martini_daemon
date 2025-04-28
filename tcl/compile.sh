#!/bin/sh

gcc -shared -o libdaemon.so daemon_tcl.c -ltcl8.6 -fPIC -g
