#!/bin/bash
fname=$1
gnuplot << _EOT_
set fit quiet
f(x) = 6.0*D*x + c
fit [0.07:] f(x) "$fname" u 1:2 via D,c
print "Diffusion Coefficient = ",D
_EOT_
