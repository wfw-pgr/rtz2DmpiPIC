import datetime
import subprocess

dlist   = [ "src", "utl", \
            "pyt", "idl"  \
        ]
name    = "stanierTube"
name    = name + str( datetime.date.today() ).replace( "-", "" )

com     = "tar zcvf " + name + ".tar.gz -X utl/tar_ex.conf"
for d in dlist:
    com =  com + " " + d
print( com )
subprocess.call( com.split(" ") )
