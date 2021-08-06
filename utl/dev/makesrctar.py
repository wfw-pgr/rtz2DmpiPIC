import datetime
import subprocess

name    = "rtz2DPIC_src"

name    = name + str( datetime.date.today() ).replace( "-", "" )
com     = "tar zcvf " + name + ".tar.gz src"
print( com )
subprocess.call( com.split(" ") )
