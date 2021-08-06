#!/bin/sh
mpiexec -machinefile utl/hosts -n 16 ./main > run.log 