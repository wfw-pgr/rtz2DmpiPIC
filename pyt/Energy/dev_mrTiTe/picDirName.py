import myStyle.LoadConfig      as lcf


def picDirName( job=None, jobDir=None ):
    if ( config is None ): config = lcf.LoadConfig()
    if ( job    is None ): job    = config["pic_job"]
    if ( jobDir is None ): jobDir = "{0}{1}".format( config["pic_jobDir"], job )
    if ( datDir is None ): datDir = "{0}dat/".format( jobDir )
    if ( pngDir is None ): pngDir = "{0}png/".format( jobDir )
    return( { "datDir":datDir, "pngDir":pngDir } )
