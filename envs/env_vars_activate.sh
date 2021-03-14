#!/bin/sh

CFLAGS_OLD=$CFLAGS
export CFLAGS_OLD
export CFLAGS="`gsl-config --cflags` ${CFLAGS_OLD}"
 
LDFLAGS_OLD=$LDFLAGS
export LDFLAGS_OLD
export LDFLAGS="`gsl-config --libs` ${LDFLAGS_OLD}"

R_LIBS_OLD=$R_LIBS
export R_LIBS_OLD
export R_LIBS=""
