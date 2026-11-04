#!/bin/bash
# usage: ./run_all.sh L   (L = 100 or 200)
# Runs the release-the-cap test (run_warm.m) and then the rho_max / p_max sweeps (run_rho.m) for one distance.
# Run from anywhere; paths follow the repository layout (code in the repository root, logs in logs/).
# The first run of this script was stopped by a time limit of the compute environment after most of the sweeps;
# run_rest.sh ran what was missing (rho_max = 1 for draw 3, p_max = 30 and 300).
L=$1
S=$(cd "$(dirname "$0")" && pwd)
cd "$S"
echo start $(date +%s) > logs/rel_start_L$L.txt
LV=$L CITERS=50 XITERS=100 DRAWS="1:4" BETAS="-5 0 5" octave-cli --no-gui run_warm.m > logs/warm_L$L.log 2>&1
echo A_done $(date +%s) >> logs/rel_start_L$L.txt
LV=$L GITERS=150 DRAWS="1:3" BETAS="0" RHOS="0.03 0.3 1.0" PMAX=100 octave-cli --no-gui run_rho.m > logs/rho_L$L.log 2>&1
LV=$L GITERS=150 DRAWS="1:3" BETAS="0" RHOS="0.1" PMAX=30 octave-cli --no-gui run_rho.m > logs/pmax30_L$L.log 2>&1
LV=$L GITERS=150 DRAWS="1:3" BETAS="0" RHOS="0.1" PMAX=300 octave-cli --no-gui run_rho.m > logs/pmax300_L$L.log 2>&1
echo B_done $(date +%s) >> logs/rel_start_L$L.txt
