import numpy                     as np
import myStyle.cMap2D            as clm
import myStyle.LoadConfig        as lcf
import myUtils.LoadFortranBinary as lfb
import myUtils.LoadConst         as lcn

def prt2cmp( job     =None , kstep   =0   , rank=0, \
             InpDir  =None , CnsFile =None, \
             keys    =None , GridMode=None, \
             x1Range =None , x2Range =None, \
             OutDir  =None , \
             silent  =False, config  =None  ):
    # --- [1] 引数チェック --- #
    if ( config   is None ): config   = lcf.LoadConfig()
    if ( job      is None ): job      = config["pic_job"]
    if ( keys     is None ): keys     = [ config["pic_key"], "psi" ]
    if ( InpDir   is None ): InpDir   = "../job/{0}/bin/kstep{1:08}/" .format( job, kstep )
    if ( CnsFile  is None ): CnsFile  = "../job/{0}/dat/constants.dat".format( job )
    if ( OutDir   is None ): OutDir   = "../job/{0}/pdf/"             .format( job )
    # --- [2] データ呼び出し --- #
    const   = lcn.LoadConst( InpFile=CnsFile )
    LI      = int( ( const["LI"] - 1 ) / 16 + 1 )
    LJ      = const["LJ"]
    print( LI, LJ )
    dct     = {"Bx":0, "ne":15, "ni":16 }
    InpFile = InpDir + "Field2D{0:08}_{1:06}.bin".format( kstep, rank )
    Data    = lfb.LoadFortranBinary( InpFile=InpFile, LI=LI, LJ=LJ, LK=18 )
    # --- [3] プロット --- #
    for key in keys:
        fData = Data[  dct[ key ] , :, : ]
        print( fData.shape )
        FigName = "{0}{1}_{2:08}_{3:06}.png".format( OutDir, key, kstep, rank )
        print( "[prt2cmp] {0} is under processing...".format( FigName ), end="" )
        fig     = clm.cMap2D( cMap=fData  , FigName=FigName )
        print( "\t [Done]" )


# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    OutDir = "./png/"
    print( "{0} is undir processing...".format( args["job"] ) )
    keys   = [ "ne"  ]
    for rank in range(16):
        drw    = prt2cmp( job =args["job"], OutDir=OutDir, keys=keys, rank=rank )
