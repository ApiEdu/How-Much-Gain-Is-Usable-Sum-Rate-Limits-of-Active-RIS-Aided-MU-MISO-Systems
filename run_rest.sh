#!/bin/bash
# usage: ./run_rest.sh L   (L = 100 or 200)
# Finishes experiment B after the first run of run_all.sh was stopped by a time limit:
# rho_max = 1.0 for draw 3 (logs/rho_L<L>_b.log), then the p_max = 30 and p_max = 300 runs (3 draws each), beta = 0 dB, 150 iterations.
L=$1
S=$(cd "$(dirname "$0")" && pwd)
cd "$S"
echo start_rest $(date +%s) > logs/rest_L$L.txt
LV=$L GITERS=150 DRAWS="3" BETAS="0" RHOS="1.0" PMAX=100 octave-cli --no-gui run_rho.m > logs/rho_L${L}_b.log 2>&1
LV=$L GITERS=150 DRAWS="1:3" BETAS="0" RHOS="0.1" PMAX=30 octave-cli --no-gui run_rho.m > logs/pmax30_L$L.log 2>&1
LV=$L GITERS=150 DRAWS="1:3" BETAS="0" RHOS="0.1" PMAX=300 octave-cli --no-gui run_rho.m > logs/pmax300_L$L.log 2>&1
echo B_done $(date +%s) >> logs/rest_L$L.txt
