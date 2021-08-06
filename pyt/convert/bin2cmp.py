import numpy                     as np
import myUtils.LoadFortranBinary as lfb
import myStyle.cMap2D            as clm
import myStyle.sfc3D             as sfc
import myStyle.LoadConfig        as lcf

def bin2cmp( InpFile=None, OutFile=None, config=None, size=None, sfcMode=False ):
    # --- [1] 引数チェック   --- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( InpFile is None ): InpFile = sys.exit("[ERROR]-@bin2cmp- InpFile=??? ")
    if ( OutFile is None ): OutFile = InpFile.replace( ".bin", ".png" )
    if ( size    is None ): size    = [ config["LI"], config["LJ"] ]
    # --- [2] データ呼び出し --- #
    print( "[bin2cmp] {0} is under processing...".format( InpFile ), end="" )
    Data        = lfb.LoadFortranBinary( InpFile=InpFile, rec=1, LI=size[0], LJ=size[1], order='F' )

    # --- [3] プロット --- #
    if ( sfcMode ):
        fig     = sfc.sfc3D ( sfc =Data, FigName=OutFile, config=config )
    else:
        fig     = clm.cMap2D( cMap=Data, FigName=OutFile, config=config )
        
    print( "\t [Done]")
    print( Data.shape )


# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args = rar.myRecvArgs()
    exe  = bin2cmp( InpFile=args["InpFile"], size=args["size"] )
