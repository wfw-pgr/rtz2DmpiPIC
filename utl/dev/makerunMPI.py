#!/usr/bin/python3
#-*- coding:utf-8 -*-

import subprocess

cmd  = 'grep -E "::\s+job\s+=" variables.f90 | cut -d \' -f2'
job  = subprocess.chech_output( cmd.split() )
