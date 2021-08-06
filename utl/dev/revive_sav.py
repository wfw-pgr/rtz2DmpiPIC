
import subprocess

# --- Input variables --- #
job   = 'tlr04-04'
step  = 46000
pnum  = 10

# --- Config --- #
dlist = [ 'dat'   ]
blist = [ 'beami', 'field', 'oldEB', 'param', 'ptspc' ]
plist = [ 'prtcl' ]
jobdir= '../job/' + job + '/sav/'

# --- .dat File Convert --- #
for subj in dlist:
    inp   = '{0}{1}_{2}_00000000.dat'.format(jobdir,job,subj)
    out   = '{0}{1}_{2}_{3:08d}.dat'.format(jobdir,job,subj,step)
    cmd   = 'mv ' + inp + ' ' + out
    print( cmd )
    # subprocess.call( cmd.split(' ') )

# --- .bin File Convert --- #
for subj in blist:
    inp   = '{0}{1}_{2}_00000000.bin'.format(jobdir,job,subj)
    out   = '{0}{1}_{2}_{3:08d}.bin'.format(jobdir,job,subj,step)
    cmd   = 'mv ' + inp + ' ' + out
    print( cmd )
    subprocess.call( cmd.split(' ') )

# --- Particle File Convert --- #
for subj in plist:
    for pn in range(pnum):
        inp   = '{0}{1}_{2}_00000000_{3:02d}.bin'.format(jobdir,job,subj,pn)
        out   = '{0}{1}_{2}_{3:08d}_{4:02d}.bin'.format(jobdir,job,subj,step,pn)
        cmd   = 'mv ' + inp + ' ' + out
        print( cmd )
        subprocess.call( cmd.split(' ') )
