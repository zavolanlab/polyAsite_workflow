#!/bin/bash

awk '{ if (NR%4==1) print ">"$0 ; if (NR%4==2) print }'
