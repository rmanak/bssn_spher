#!/bin/bash

# Set up your maple executable directory
MAPLE=/opt/maple18/bin/maple

VERBOSE=0

echo "using: "
echo "maple   -> $MAPLE"
echo "verbose -> $VERBOSE"
echo "==================================="

which $MAPLE > /dev/null 2>&1 || echo "cannot find $MAPLE"
which $MAPLE > /dev/null 2>&1 || exit 1



echo "deriving BSSN, performing compactification and FDA, generating solver and IRE codes..."

if [ $VERBOSE = 1 ]; then
	$MAPLE < bssn_main.mpl
else
	$MAPLE < bssn_main.mpl | tail -1
fi
echo "done!"
echo "=================================="

echo "building FDA codes for coordinate choices, ctfm initializer, apparent horizon finder "
if [ $VERBOSE = 1 ]; then
	$MAPLE < coord_choices.mpl
	$MAPLE < ctfm.mpl
	$MAPLE < ctfm2.mpl
	$MAPLE < trap_surface.mpl
else
	$MAPLE < coord_choices.mpl | tail -1
	$MAPLE < ctfm.mpl | tail -1
	$MAPLE < ctfm2.mpl | tail -1
	$MAPLE < trap_surface.mpl | tail -1
fi
echo "done!"
echo "=================================="

/bin/mv *.f *.h *_call ./src/FD/
