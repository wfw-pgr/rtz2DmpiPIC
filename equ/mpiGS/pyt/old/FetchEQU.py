import numpy                     as np
import myUtils.LoadFortranBinary as lfb
import myUtils.LoadConst         as lcn
import myStyle.LoadConfig        as lcf
import myConvert.pilearr         as pil
import PICsees.calcFlux          as cfl
import fLIB.fLIB__rct2rct        as r2r

# -- Flux は [LJ,LI], Data は [LI,LJ] で返却されることに注意 -- #

def FetchEQU( job     =None , \
              InpDir  =None , CnsFile =None, \
              keys    =None , GridMode=None, \
              x1Range =None , x2Range =None, \
              silent  =False, config  =None \
):
    # --- [1] 引数チェック --- #
    if ( config   is None ): config   = lcf.LoadConfig()
    if ( job      is None ): job      = config["equ_job"]
    if ( GridMode is None ): GridMode = config["equ_GridMode"]
    if ( keys     is None ): keys     = [ config["equ_key"], "psi" ]
    if ( InpDir   is None ): InpDir   = "../job/{0}/{1}/"             .format( job, GridMode )
    if ( CnsFile  is None ): CnsFile  = "../job/{0}/{1}/parameter.dat".format( job, GridMode )
    
    # --- [2] データ呼び出し --- #
    const   = lcn.LoadConst( InpFile=CnsFile )
    xAxis   = np.linspace( const["x1Min"], const["x1Max"], const["LI"] )
    yAxis   = np.linspace( const["x2Min"], const["x2Max"], const["LJ"] )
    #  -- [2-1] Bfd -- #
    InpFile = InpDir + "Bfd.bin"
    Bfd     = lfb.LoadFortranBinary( InpFile =InpFile,     \
                                     rec     =1,           \
                                     LI      =3,           \
                                     LJ      =const["LI"], \
                                     LK      =const["LJ"], \
                                     order   =config["LoadBinaryOrder"] \
                                 )
    #  -- [2-2] Jcr -- #
    InpFile = InpDir + "Jcr.bin"
    Jcr     = lfb.LoadFortranBinary( InpFile =InpFile,     \
                                     rec     =1,           \
                                     LI      =3,           \
                                     LJ      =const["LI"], \
                                     LK      =const["LJ"], \
                                     order   =config["LoadBinaryOrder"] \
                                 )
    #  -- [2-3] frp -- #
    InpFile = InpDir + "frp.bin"
    frp     = lfb.LoadFortranBinary( InpFile =InpFile,     \
                                     rec     =1,           \
                                     LI      =3,           \
                                     LJ      =const["LI"], \
                                     LK      =const["LJ"], \
                                     order   =config["LoadBinaryOrder"] \
                                 )
    #  -- [2-4] uvc -- #
    InpFile = InpDir + "uvc.bin"
    uvc     = lfb.LoadFortranBinary( InpFile =InpFile,     \
                                     rec     =1,           \
                                     LI      =3,           \
                                     LJ      =const["LI"], \
                                     LK      =const["LJ"], \
                                     order   =config["LoadBinaryOrder"] \
                                 )
    Data    = pil.pilearr( ( Bfd.transpose(2,0,1), Jcr.transpose(2,0,1),
                             frp.transpose(2,0,1), uvc.transpose(2,0,1) ) \
                           , NoNewAxis=True )
    # --- [3] Flux 磁束関数の返却 --- #
    flux = cfl.calcFlux( rg=yAxis, Bz=Data[2,:,:], indexing='ji' )
    if ( config["reverseFlux"]   ): flux = - flux
    if ( config["NormalizeFlux"] ):
        LIc, LJc = int( flux.shape[1]/2), int( flux.shape[0]/2 )
        if ( flux[LJc,LIc] > 0.0 ): flux = flux / np.max( flux )
        if ( flux[LJc,LIc] < 0.0 ): flux = flux / np.min( flux )

    # --- [4] 切り出し : xRange, yRange で 検索して 返却 -- #
    if ( x1Range is not None ):
        idx   = ( np.where( ( xAxis > x1Range[0] ) & ( xAxis < x1Range[1] ) ) )[0]
        Data  =  Data[:,:,idx]
        xAxis = xAxis[    idx]
        if ( config["Data_retPsi"] ): flux = flux[:,idx]
    if ( x2Range is not None ):
        idx   = ( np.where( ( yAxis > x2Range[0] ) & ( yAxis < x2Range[1] ) ) )[0]
        Data  =  Data[:,idx,:]
        yAxis = yAxis[    idx]
        if ( config["Data_retPsi"] ): flux = flux[idx,:]

    # --- [5] 返却準備 --- #
    #  -- [5-1] 対応表 --  #
    tbl   = { "Bx" :0, "By" :1, "Bz" :2, "Jx":3, "Jy":4 , "Jz":5,
              "rho":7, "prs":8, "ux":9, "uy":10, "uz":11 }
    #  -- [5-2] 荷造り -- #
    pack  = '{ "xAxis":xAxis, "yAxis":yAxis, "flux":flux'
    for key in keys:
        if  ( key in tbl.keys() ):
            pack  += ', "' + str(key) + '":Data[{0},:,:]'.format( tbl[key] )
        elif( key == "psi"      ):
            pack  += ', "psi":flux'
        else:
            sys.exit( " [ERROR] Unknown key -@FetchEQU- [ERROR] " )
    pack  = pack + "}"
    
    # --- [6] 返却 --- #
    if ( not(silent) ):
        print()
        print( "[FetchEQU] :: Return" )
        print( "\t {0}".format(pack) )
    ret   = eval( pack )
    return( ret )



# ----------------------------------------- #
if ( __name__=="__main__" ):
    # -- 引数 -- #
    import myUtils.myRecvArgs as rar
    args      = rar.myRecvArgs()
    Fields  = FetchEQU( keys=["Bx", "By", "Bz"] )
    print( Fields["Bx"].shape ) 


