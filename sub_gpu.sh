#!/bin/bash
#SBATCH --export=ALL
#SBATCH --job-name=298_0.1
#SBATCH -o out -e err
#SBATCH -N 1 -n 1 -t 300:00:00 --mem 10000mb -c 1
#SBATCH --gres=gpu:1 -p gpu

./run_gcmc.py > logfile

