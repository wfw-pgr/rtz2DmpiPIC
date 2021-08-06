import datetime
import subprocess

dlist   = [ "src", "utl",
            "equ", "pyt"
]
name    = "rtz2DmpiPIC"

name    = name + str( datetime.date.today() ).replace( "-", "" )
com     = "tar zcvf " + name + ".tar.gz -X utl/tar_ex.conf"
for d in dlist:
    com =  com + " " + d
print( com )
subprocess.call( com.split(" ") )
