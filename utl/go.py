import os, sys, subprocess, argparse

# --- 引数 --- #
parser = argparse.ArgumentParser()
parser.add_argument("-nolog", help="the init time step.",action="store_true")
args   = parser.parse_args()

# --- 設定 --- #
HostsF = "utl/hosts"
nProcs = 0

# --- Host 情報の取得 --- #
with open( HostsF, 'r' ) as hf:
    for line in hf.readlines():
        ppp    = int( ( line.split(':')[1] ).strip() )
        nProcs = nProcs + ppp                  # - 総プロセッサ数の取得 -- #
nHosts   = sum( 1 for i in open(HostsF, 'r') ) # - HostsFの行数取得 - # 

# --- Shell ファイル作成 及び 実行 --- #
shellhdr = "#!/bin/sh\n"
cmd      = "mpiexec -machinefile {0} -n {1} ./main".format( HostsF, nProcs )
if not( args.nolog ):
    cmd  = shellhdr + cmd + " > run.log "
    with open( "utl/run.sh", "w" ) as rh:
        rh.write( cmd )
    exe  = subprocess.check_output( "sh utl/run.sh".split() )
else:
    exe  = subprocess.call( cmd.split() )
