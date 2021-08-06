import sys
import numpy                     as np
import PICsees.FetchPIC          as fpc
import PICsees.createMask        as msk
import myStyle.LoadConfig        as lcf
import myStyle.cMap2D            as clm
import myUtils.LoadConst         as lcn
import myConvert.pilearr         as pil
import fLIB.fLIB__LinearFilter2D as flt
import fLIB.fLIB__curl           as crl
import fLIB.fLIB__ftranspose     as ftr
import fLIB.fLIB__helicity       as ghl
import myConvert.where2D         as wh2
import myStyle.plot1D            as pl1
import myBasicAlgs.getPsiMap     as gPM

# --- ヘリシティを算出する --- #
def pic2hlx( InpFile=None, InpDir=None, OutFile=None , OutDir=None, \
             CnsFile=None, config=None ):
    # --- [1] 引数チェック --- #
    if ( config  is     None ): config   = lcf.LoadConfig()
    if ( CnsFile is     None ): CnsFile  = "../job/{0}/dat/constants.dat".format( config["job"] )
    if ( InpFile is     None ): InpFile  = "Field2D{0:08}.bin".format( config["kstep"] )
    if ( InpFile is     None ): InpDir   =   "../job/{0}/bin/".format( config["job"  ] )
    if ( OutFile is     None ): OutFile  = InpFile.replace( ".bin", ".png" )
    if ( OutDir  is not None ): OutFile  = OutDir + OutFile
    # --- [2] データ呼び出し --- #
    #  -- [2-1] Data Fetch  --  #
    const       = lcn.LoadConst( InpFile=CnsFile )
    keys        = ["Bx","By","Bz","uix","uiy","uiz","psi"]
    Data        = fpc.FetchPIC( InpFile=InpFile, InpDir=InpDir, \
                                CnsFile=CnsFile, keys  =keys  , \
                                config =config )
    zAxis       = Data["xAxis"]
    rAxis       = Data["yAxis"]
    flux        = Data["psi"  ]
    uvc         = pil.pilearr( ( ftr.ftranspose( Data["uix"] ), \
                                 ftr.ftranspose( Data["uiz"] ), \
                                 ftr.ftranspose( Data["uiy"] )  ) )
    bfd         = pil.pilearr( ( ftr.ftranspose( Data["Bx" ] ), \
                                 ftr.ftranspose( Data["Bz" ] ), \
                                 ftr.ftranspose( Data["By" ] )  ) )
    wvc         = crl.curlv  ( u=uvc, x1=rAxis, x2=zAxis, coordinate="rzt" )
    avp         = crl.curlInv( b=bfd, x1=rAxis, x2=zAxis, coordinate="rzt" )
    volume      = const["dr"] * const["dz"]

    # -- hlx -- #
    #  - 1 - #
    hlxSimTot   = ghl.GridHelicity( u=uvc, w=wvc, a=avp, b=bfd, vol=volume )
    #  - 2 - #
    Mask        = np.transpose( msk.createMask( flux=flux, criterion=0.05 ) )
    MaskSep     = pil.pilearr( (Mask,Mask,Mask) )
    uvcInSep    = uvc * MaskSep
    wvcInSep    = wvc * MaskSep
    bfdInSep    = bfd * MaskSep
    avpInSep    = avp * MaskSep
    hlxInSep    = ghl.GridHelicity( u=uvcInSep, w=wvcInSep, a=avpInSep, b=bfdInSep, vol=volume )
    #  - 3 - #
    xptRange    = [0.4,0.2,0.6,0.8]
    psiXO       = gPM.getPsiMap( psi     =flux, x1g=zAxis, x2g=rAxis, \
                                 xptRange=xptRange )
    privMask    = np.transpose( msk.createMask( flux=flux, criterion=psiXO["psiX"] ) )
    MaskUps     = ( wh2.where2D( Data   =privMask, x1Axis=zAxis, x1Range=[psiXO["x1Xpt"],zAxis[-1]], \
                                 Masking=True ) )["Data"]
    MaskUp      = pil.pilearr( (MaskUps,MaskUps,MaskUps) )
    uvcUp       = uvc * MaskUp
    wvcUp       = wvc * MaskUp
    bfdUp       = bfd * MaskUp
    avpUp       = avp * MaskUp
    hlxUp       = ghl.GridHelicity( u=uvcUp   , w=wvcUp   , a=avpUp   , b=bfdUp   , vol=volume )
    #  - 4 - #
    MaskDws     = ( wh2.where2D( Data   =privMask, x1Axis=zAxis, x1Range=[zAxis[ 0],psiXO["x1Xpt"]], \
                                 Masking=True ) )["Data"]
    MaskDw      = pil.pilearr( (MaskDws,MaskDws,MaskDws) )
    uvcDw       = uvc * MaskDw
    wvcDw       = wvc * MaskDw
    bfdDw       = bfd * MaskDw
    avpDw       = avp * MaskDw
    hlxDw       = ghl.GridHelicity( u=uvcDw   , w=wvcDw   , a=avpDw   , b=bfdDw   , vol=volume )    
    # --- 返却 --- #
    ret         = [ hlxSimTot, hlxInSep, hlxUp, hlxDw ]
    return( ret )


def pic2hlx_Driver( InpDir=None, ksteps=None, OutFile=None, OutDir=None, CnsFile=None, config=None ):
    # -- 引数チェック -- #
    if ( config  is     None ): config  = lcf.LoadConfig()
    if ( CnsFile is     None ): CnsFile  = "../job/{0}/dat/constants.dat".format( config["job"] )
    if ( OutFile is     None ): OutFile  = "picHelicity.dat"
    if ( OutDir  is not None ): OutFile  = OutDir + OutFile
    if ( ksteps  is     None ): ksteps   = [0]
    f = open( OutFile,"w" )
    f.close()
    
    for kstep in ksteps:
        binFile = "Field2D" + "{0:08d}".format( kstep ) + ".bin"
        print( "[pic2hlx] {0} is processing...".format( binFile ), end="" )
        ret     = pic2hlx( InpFile=binFile, InpDir=InpDir, OutDir=OutDir,
                           CnsFile=CnsFile, config=config )
        print( "\t [completed]" )
        HelixV  = np.concatenate( [ ( ret[0] )["HelixV"], ( ret[1] )["HelixV"],
                                    ( ret[2] )["HelixV"], ( ret[3] )["HelixV"] ] )
        with open( OutFile, "ab" ) as f:
            np.savetxt( f, HelixV[np.newaxis,:] )
        print( HelixV[:] )

def dat2fig( datFile=None, datDir=None, OutDir=None, OutFile=None, CnsFile=None, config=None ):
    # -- 引数チェック -- #
    if ( config  is     None ): config  = lcf.LoadConfig()
    if ( CnsFile is     None ): CnsFile  = "../job/{0}/dat/constants.dat".format( config["job"] )
    if ( datFile is     None ): datFile  = "picHelicity.dat"
    if ( OutFile is     None ): OutFile  = "picHelicity.png"
    if ( datDir  is not None ): datFile  = datDir + datFile
    if ( OutDir  is not None ): OutFile  = OutDir + OutFile
    const       = lcn.LoadConst( InpFile=CnsFile )
    Data        = np.transpose( np.loadtxt( datFile, comments="#" ) )
    tAxis       = const["dt"] * np.arange( Data.shape[1] ) * 1000.

    # --- [2] コンフィグ --- #
    config["plt_LegNColumn"] = 3
    config["FigSize"]        = (6,2)
    config["xTitle"]         = "$t$ / $t_A$"
    config["plt_AutoRange"]  = False
    config["AutoTicks"]      = False
    config["plt_xRange"]     = [0,4000]
    config["plt_yRange"]     = [-5,5]
    config["plt_SplineItp"]  = True
    config["plt_LegFontSize"]= 10
    config["plt_Position"]   = [0.16,0.25,0.97,0.97]

    # --- [2] シミュレーション領域 --- #
    config["yTitle"]         = ""
    config["plt_AutoRange"]  = True
    config["AutoRange"]      = True
    config["AutoTicks"]      = True
    OutFile                  = OutDir + "hlxInSim.png"
    img                      = pl1.plot1D( FigName=OutFile, config=config )
    y1 = Data[0,:] * const["rmi"]**2
    y2 = Data[1,:] * const["rmi"] *  2.0
    y3 = Data[2,:]
    img.addPlot( xAxis=tAxis, yAxis=y1, label="Fluid" )
    img.addPlot( xAxis=tAxis, yAxis=y2, label="Cross" )
    img.addPlot( xAxis=tAxis, yAxis=y3, label="Magnetic" )
    img.addLegend()
    img.writeFigure()
    print( "[pic2hlx] Figure : {0} is outputed... ".format(OutFile) )

    # --- [2] セパラトリクス領域 --- #
    config["yTitle"]         = ""
    config["plt_AutoRange"]  = False
    config["AutoRange"]      = False
    config["AutoTicks"]      = True
    OutFile                  = OutDir + "hlxInSep.png"
    img                      = pl1.plot1D( FigName=OutFile, config=config )
    y1 = Data[3,:] * const["rmi"]**2
    y2 = Data[4,:] * const["rmi"] *  2.0
    y3 = Data[5,:]
    print( const["rmi"] )
    img.addPlot( xAxis=tAxis, yAxis=y1, label="Fluid" )
    img.addPlot( xAxis=tAxis, yAxis=y2, label="Cross" )
    img.addPlot( xAxis=tAxis, yAxis=y3, label="Magnetic" )
    img.addLegend()
    img.writeFigure()
    print( "[pic2hlx] Figure : {0} is outputed... ".format(OutFile) )

    # --- [3] プライベート領域 --- #
    config["yTitle"]         = ""
    config["plt_AutoRange"]  = False
    config["AutoRange"]      = False
    config["plt_xRange"]     = [0,4000]
    config["plt_yRange"]     = [-5,5]
    config["AutoTicks"]      = True
    OutFile                  = OutDir + "hlxPrivate.png"
    img                      = pl1.plot1D( FigName=OutFile, config=config )
    y1 = Data[ 6,:] * const["rmi"]**2
    y2 = Data[ 7,:] * const["rmi"] *  2.0
    y3 = Data[ 8,:]
    y4 = Data[ 9,:] * const["rmi"]**2
    y5 = Data[10,:] * const["rmi"] *  2.0
    y6 = Data[11,:]
    print( const["rmi"] )
    img.addPlot( xAxis=tAxis, yAxis=y1, label="Fluid #1" )
    img.addPlot( xAxis=tAxis, yAxis=y2, label="Cross #1" )
    img.addPlot( xAxis=tAxis, yAxis=y3, label="Magnetic #1" )
    img.addPlot( xAxis=tAxis, yAxis=y4, label="Fluid #2" )
    img.addPlot( xAxis=tAxis, yAxis=y5, label="Cross #2" )
    img.addPlot( xAxis=tAxis, yAxis=y6, label="Magnetic #2" )
    img.addLegend()
    img.writeFigure()
    print( "[pic2hlx] Figure : {0} is outputed... ".format(OutFile) )

    
        
# --------------------------------------- #
if ( __name__=="__main__" ):
    # --  引数   -- #
    import myUtils.myRecvArgs as rar
    args        = rar.myRecvArgs()
    CnsFile     = args["jobDir"] + "dat/constants.dat"
    InpDir      = args["jobDir"] + "bin/"
    OutDir      = args["jobDir"] + "png/"
    OutDir      = "./png/"
    config      = lcf.LoadConfig()
    datDir      = "./dat/"
    datFile     = "picHelicity.dat"

    # --  解析   -- #
    Flag__pic2hlxAnalyzer = False
    if ( Flag__pic2hlxAnalyzer ):
        exe = pic2hlx_Driver( InpDir =InpDir , ksteps=args["ksteps"], OutFile=datFile, OutDir=datDir, \
                              CnsFile=CnsFile, config=config )

    # -- お絵かき -- #
    Flag__Figure = True
    if ( Flag__Figure ):
        exe = dat2fig( CnsFile=CnsFile, datDir=datDir, datFile=datFile, OutDir=OutDir, config=config )
