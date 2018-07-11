#!/bin/sh

# key                 v1in   v2in  v3in  v1odd  nevents ewights pweights conservepT
root -l -q -b 'v1MC+("0.000, 0.06, 0.03, false, 1e6,    true,   false,   false")'
root -l -q -b 'v1MC+("0.000, 0.06, 0.03, false, 1e6,    true,   false,   true")'
root -l -q -b 'v1MC+("0.000, 0.06, 0.03, false, 1e6,    true,   true,    true")'

root -l -q -b 'v1MC+("0.015, 0.06, 0.03, false, 1e6,    true,   false,   false")'
root -l -q -b 'v1MC+("0.015, 0.06, 0.03, false, 1e6,    true,   true,    false")'
root -l -q -b 'v1MC+("0.015, 0.06, 0.03, false, 1e6,    true,   false,   true")'
root -l -q -b 'v1MC+("0.015, 0.06, 0.03, false, 1e6,    true,   true,    true")'

root -l -q -b 'v1MC+("0.015, 0.06, 0.03, true,  1e6,    true,   false,   false")'
root -l -q -b 'v1MC+("0.015, 0.06, 0.03, true,  1e6,    true,   false,   true")'
