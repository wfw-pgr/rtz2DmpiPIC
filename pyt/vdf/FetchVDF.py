import numpy                        as np
import myStyle.LoadConfig           as lcf
import IORoutines.LoadFortranBinary as lfb

# ================================================================ #
# ===  FetchVDF :: VDFをロードして返却                         === #
# ================================================================ #
def FetchVDF( job   =None, kstep  =None, jobDir =None, datDir=None, config=None, returnType="direct", \
              binDir=None, vdfFile=None, vcnFile=None, vdfset=None, axtype=None, nstype=None ):
    # ---------------------------------------- #
    # --- [1]    引数チェック              --- #
    # ---------------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( datDir  is None ): datDir  = "{0}dat/".format( jobDir )
    if ( binDir  is None ): binDir  = "{0}bin/kstep{1:08}/"    .format( jobDir, kstep )
    if ( vdfFile is None ): vdfFile = "{0}VDFsample_{1:08}.bin".format( binDir, kstep )
    if ( vcnFile is None ): vcnFile = "{0}vdfConst.dat".format( datDir )
    if ( vdfset  is None ): vdfset  = 0
    if ( axtype  is None ): axtype  = "xy"
    if ( nstype  is None ): nstype  = "i"
    # ---------------------------------------- #
    # --- [2]    読込                      --- #
    # ---------------------------------------- #
    with open( vcnFile, "r" ) as f:
        line1 = ( f.readline() ).split()
        line2 = ( f.readline() ).split()
        Ndiv  = int( line1[1] )
        Nset  = int( line2[1] )
    vdfcns    = ( np.loadtxt( vcnFile, skiprows=3 ) )
    size      = ( Nset, 2, 3, Ndiv+2, Ndiv+2 )
    vdf       = lfb.LoadFortranBinary( FileName=vdfFile, size=size )
    # ---------------------------------------- #
    # --- [3]    返却                      --- #
    # ---------------------------------------- #
    if   ( returnType == "direct" ):
        ret   = { "Ndiv":Ndiv, "Nset":Nset, "vdfcns":vdfcns, "vdf":vdf }
    elif ( returnType == "select" ):
        ax    = { "xy":0, "yz":1, "zx":2 }[axtype]
        ns    = { "e" :0, "i" :1         }[nstype]
        ret   = { "Ndiv":Ndiv, "Nset":Nset, "vdfcns":vdfcns, "vdf":vdf[:,:,ax,ns,vdfset] }
    return( ret )


# ================================================================ #
# ===  実行部                                                  === #
# ================================================================ #
if ( __name__=="__main__" ):
    import myUtils.genArgs as gar
    args    = gar.genArgs()
    vdf     = FetchVDF( job=args["job"], jobDir=args["jobDir"] )
    print( vdf["Ndiv"], vdf["Nset"], vdf["vdfcns"].shape, vdf["vdf"].shape )
