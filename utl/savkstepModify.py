import subprocess
import myUtils.getFileList as gfl

def savkstepModify( Execute=True, manual=False ):
    # --------------- #
    # -- parameter -- #
    # --------------- #
    if ( manual ):
        job   = "CtrO22_01"
        bfr   = 0
        aft   = 284000
        PEtot = 16
    else:
        File  = gfl.getFileList( search="const/*jobcg_00000000.dat" )
        with open( File[0], "r" ) as f:
            lines = f.readlines()
        job   = str( lines[0] ).strip()
        bfr   = 0
        aft   = int( lines[2] )
        PEtot = 16
    # --------------- #
    # --   prtcl   -- #
    # --------------- #
    tgt   = "prtcl"
    for iPE in range( PEtot ):
        f1  = "{1}/{0}_{1}_{2:08}_{3:06}.bin".format( job, tgt, bfr, iPE )
        f2  = "{1}/{0}_{1}_{2:08}_{3:06}.bin".format( job, tgt, aft, iPE )
        odr = "mv {0} {1}".format( f1, f2 )
        print( odr )
        if ( Execute ): subprocess.call( odr.split() )
    # --------------- #
    # --   field   -- #
    # --------------- #
    tgt   = "field"
    for iPE in range( PEtot ):
        f1  = "{1}/{0}_{1}_{2:08}_{3:06}.bin".format( job, tgt, bfr, iPE )
        f2  = "{1}/{0}_{1}_{2:08}_{3:06}.bin".format( job, tgt, aft, iPE )
        odr = "mv {0} {1}".format( f1, f2 )
        print( odr )
        if ( Execute ): subprocess.call( odr.split() )
    # --------------- #
    # --   oldEB   -- #
    # --------------- #
    tgt   = "oldEB"
    for iPE in range( PEtot ):
        f1  = "{1}/{0}_{1}_{2:08}_{3:06}.bin".format( job, tgt, bfr, iPE )
        f2  = "{1}/{0}_{1}_{2:08}_{3:06}.bin".format( job, tgt, aft, iPE )
        odr = "mv {0} {1}".format( f1, f2 )
        print( odr )
        if ( Execute ): subprocess.call( odr.split() )
    # --------------- #
    # --   mpicg   -- #
    # --------------- #
    tgt   = "mpicg"
    for iPE in range( PEtot ):
        f1  = "{1}/{0}_{1}_{2:08}_{3:06}.dat".format( job, tgt, bfr, iPE )
        f2  = "{1}/{0}_{1}_{2:08}_{3:06}.dat".format( job, tgt, aft, iPE )
        odr = "mv {0} {1}".format( f1, f2 )
        print( odr )
        if ( Execute ): subprocess.call( odr.split() )
    # --------------- #
    # --   const   -- #
    # --------------- #
    #  - jobcg - #
    tgt = "jobcg"
    f1  = "const/{0}_{1}_{2:08}.dat".format( job, tgt, bfr )
    f2  = "const/{0}_{1}_{2:08}.dat".format( job, tgt, aft )
    odr = "mv {0} {1}".format( f1, f2 )
    print( odr )
    if ( Execute ): subprocess.call( odr.split() )
    #  - param - #
    tgt = "param"
    f1  = "const/{0}_{1}_{2:08}.bin".format( job, tgt, bfr )
    f2  = "const/{0}_{1}_{2:08}.bin".format( job, tgt, aft )
    odr = "mv {0} {1}".format( f1, f2 )
    print( odr )
    if ( Execute ): subprocess.call( odr.split() )
    #  - ptspc - #
    tgt = "ptspc"
    f1  = "const/{0}_{1}_{2:08}.bin".format( job, tgt, bfr )
    f2  = "const/{0}_{1}_{2:08}.bin".format( job, tgt, aft )
    odr = "mv {0} {1}".format( f1, f2 )
    print( odr )
    if ( Execute ): subprocess.call( odr.split() )
    #  - ijDomain - #
    f1  = "const/ijDomain{0:08}.dat".format( bfr )
    f2  = "const/ijDomain{0:08}.dat".format( aft )
    odr = "mv {0} {1}".format( f1, f2 )
    print( odr )
    if ( Execute ): subprocess.call( odr.split() )


if ( __name__=="__main__" ):
    savkstepModify()
