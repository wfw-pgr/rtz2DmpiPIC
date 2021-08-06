import re, sys
import numpy               as np
import myBasicAlgs.getPsiX as gPX
import myStyle.myMPL       as mml
import myStyle.LoadConfig  as lcf
import myUtils.LoadConst   as lcn
import PICsees.FetchField  as ffl
import myConvert.where2D   as wh2
import myConvert.argmax2D  as am2
import myConvert.pilearr   as pil

# --- リコネクションレートの計算材料 を 計算 --- #
def RecRateMaterial( InpFile=None,        InpDir =None,        CnsFile=None , \
                     zXRange=[0.40,0.60], rXRange=[0.20,0.80], silent =False, \
                     config =None ):
    # -- [0] 引数チェック -- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( CnsFile is None ): CnsFile = "../job/{0}/dat/constants.dat".format( config["job"] )
    if ( InpFile is None ):
        InpFile =  "Field2D{0:08}.bin".format( config["kstep"] )
        InpDir  = "../../job/{0}/bin/".format( config["job"  ] )

    # -- [1] データ読込 -- #
    keys        = ["Bx", "By", "ni"]
    const       = lcn.LoadConst ( InpFile=CnsFile  )
    Field       = ffl.FetchField( InpFile=InpFile, InpDir =InpDir , \
                                  CnsFile=CnsFile, keys   =keys   , silent =True )
    zAxis       = Field["xAxis"]
    rAxis       = Field["yAxis"]
    flux        = Field["flux" ]
    kstep       = int( ( re.search( "[0-9]{8}", InpFile ) ).group(0) )
    time        = const["dt"] * float( kstep )

    # -- [2] X点での psi -- #
    config["NormalizeFlux"] = False
    psiX,iX,jX  = gPX.getPsiX( psi=flux, x1Range=rXRange, x2Range=zXRange )
    #  - 鞍点みつからないとき - #
    if ( psiX is None ):
        psiX    = np.max( flux ) # - 最大値 - #
        jX, iX  = am2.argmax2D( flux )
    zXpt, rXpt  = zAxis[iX], rAxis[jX]

    # -- [3] O点での psi -- #
    #  - [3-1] 左 - #
    psiF1d      = wh2.where2D( Data=flux, x2Range=[0.0,0.5] )
    psiO1       = np.max( psiF1d["Data"] ) # - 最大値 - #
    ijO1        = am2.argmax2D( psiF1d["Data"] )
    jO1         = ( psiF1d["LBRT"] )[0]+ijO1[0]
    iO1         = ( psiF1d["LBRT"] )[1]+ijO1[1]
    rOpt1,zOpt1 = rAxis[jO1], zAxis[iO1]
    #  - [3-2] 右 - #
    psiF2d      = wh2.where2D( Data=flux, x2Range=[0.5,1.0] )
    psiO2       = np.max( psiF2d["Data"] ) # - 最大値 - #
    ijO2        = am2.argmax2D( psiF2d["Data"] )
    jO2         = ( psiF2d["LBRT"] )[0]+ijO2[0]
    iO2         = ( psiF2d["LBRT"] )[1]+ijO2[1]
    rOpt2,zOpt2 = rAxis[jO2], zAxis[iO2]

    # -- [4]   レートの算出   -- #
    #  - [4-1] O点間で最強のb_m, 及び，そこでのrho_mを取得 - #
    jO         =      int( 0.5*( jO1+jO2 ) )
    bx_        = ( Field["Bx"] )[jO,iO1:iO2]
    by_        = ( Field["By"] )[jO,iO1:iO2]
    ni_        = ( Field["ni"] )[jO,iO1:iO2]
    absB       = np.sqrt( bx_**2 + by_**2 )
    im         = np.argmax( absB )
    bm         =              absB[ im ]
    rhom       = const["mr"] * ni_[ im ]
    vAm        = bm / ( const["wpewce"] * np.sqrt( rhom ) )
    # --- [5] レートの計算 材料 書き出し --- #
    ret        = np.array( [kstep, time, bm, vAm, rXpt, psiO1, psiO2, psiX] )
    return( ret )


# --- リコネクションレート を計算 --- #
def RecRateMaterialDriver( ksteps =None     , InpDir ="./", \
                           OutDir =None     , OutFile='RecRateMaterials.dat', \
                           x1Range=[0.4,0.6], x2Range=[0.2,0.8], \
                           CnsFile=None ):
    if ( ksteps is     None ): sys.exit( "[ERROR] No ksteps -@RecRateMateDriver- [ERROR]" )
    if ( OutDir is not None ): OutFile = OutDir + OutFile
    for kstep in ksteps:
        InpFile   = "Field2D{0:08}.bin".format( kstep )
        print( " RecRate :: {0} is under processing...".format( InpFile ) )
        RecMate = RecRateMaterial( InpFile=InpFile, \
                                   InpDir =InpDir , \
                                   CnsFile=CnsFile, \
                                   zXRange=x1Range, \
                                   rXRange=x2Range  )
        with open( OutFile, "ab" ) as f:
            np.savetxt( f, RecMate[np.newaxis,:] )

        
def RecRate( InpFile='RecRateMaterials.dat', OutFile='RecRate.dat' ):
    Mats   = np.loadtxt( InpFile )
    print( Mats.shape )
    Nfiles = Mats.shape[0]
    time   = np.zeros( (Nfiles,) )
    absB   = np.zeros( (Nfiles,) )
    vAin   = np.zeros( (Nfiles,) )
    dpsi   = np.zeros( (Nfiles,) )
    rXpt   = np.zeros( (Nfiles,) )
    dpdt   = np.zeros( (Nfiles,) )
    RecR   = np.zeros( (Nfiles,) )
    psiX   = np.zeros( (Nfiles,) )
    mrgR   = np.zeros( (Nfiles,) )
    for i in range( Nfiles ):
        psiOpt  = np.max( [Mats[i,5], Mats[i,6]] )
        dpsi[i] = psiOpt - Mats[i,7]
        psiX[i] = Mats[i,7] / psiOpt
        if ( abs( dpsi[i] ) < abs( 0.05*psiOpt ) ):
            dpsi[i] = 0.0
    for i in range( 1,Nfiles-1 ):
        dtInv   = 1.0/( Mats[i+1,1]-Mats[i,1] )
        time[i] = 0.5*( Mats[i+1,1]+Mats[i,1] )
        absB[i] = 0.5*( Mats[i+1,2]+Mats[i,2] )
        vAin[i] = 0.5*( Mats[i+1,3]+Mats[i,3] )
        rXpt[i] = 0.5*( Mats[i+1,4]+Mats[i,4] )
        mrgR[i] = 0.5*( psiX[i+1  ]+psiX[i  ] )
        dpdt[i] = np.abs( ( dpsi[i+1]-dpsi[i] )*dtInv )
        RecR[i] = dpdt[i] / ( rXpt[i]*absB[i]*vAin[i] )
    Rates  = pil.pilearr( (time,absB,vAin,rXpt,dpsi,RecR,mrgR) )
    np.savetxt( OutFile, Rates, delimiter=' ' )

    
# --- Reconection Rate をプロット --- #
def RecRatePlotter( InpFile='RecRate.dat', OutDir="./png/", config=None ):
    # -- 設定情報 / 読込 -- #
    if ( config is None ): config = lcf.LoadConfig()
    data   = np.loadtxt( InpFile )
    ndata  = len( data[0,:] )
    # -- [1] プロットのコンフィグ -- #
    config["AutoRange"] = False
    # -- [2] リコネクション率 -- #
    xAxis  = data[0,1:-1]
    yAxis  = data[5,1:-1]
    xRange = [0.0, 2000.0]
    yRange = [+0.0,+2.e-1]
    OutFile= OutDir + "RecRate.png"
    fig    = mml.plot1D( xAxis  =xAxis  , yAxis  =yAxis   , \
                         xRange =xRange , yRange =yRange  , \
                         xTitle ='Time' , yTitle =r'$E_R$', \
                         FigSize=(5,2.5), FigName=OutFile , \
                         config =config                       )
    # -- [3] X-point R 座標 -- #
    xAxis  = data[0,1:-1]
    yAxis  = data[3,1:-1]
    xRange = [0.0, 2000.0]
    yRange = [+0.0,+2.e+1]
    OutFile= OutDir + "rXpoint.png"
    fig    = mml.plot1D( xAxis  =xAxis  , yAxis  =yAxis   , \
                         xRange =xRange , yRange =yRange  , \
                         xTitle ='Time' , yTitle =r'$R_X$', \
                         FigSize=(5,2.5), FigName=OutFile , \
                         config =config                       )
    # -- [4] 合体率 -- #
    xAxis  = data[0,1:-1]
    yAxis  = data[6,1:-1] * 100.0
    xRange = [0.0, 2000.0]
    yRange = [+0.0,+120.0]
    OutFile= OutDir + "MergingRate.png"
    fig    = mml.plot1D( xAxis  =xAxis  , yAxis  =yAxis   , \
                         xRange =xRange , yRange =yRange  , \
                         xTitle ='Time' , yTitle =r'$\eta_{M}$', \
                         FigSize=(5,2.5), FigName=OutFile , \
                         config =config                       )
    # -- [5] 合体率 / リコネクション率 を両方描く -- #
    xAxis      = data[0,1:]
    yAxis1     = data[5,1:]
    yAxis2     = data[6,1:] * 0.15
    xAxis [-1] = 3000.
    yAxis1[-1] = yAxis1[-2]
    yAxis2[-1] = yAxis2[-2]
    xRange = [0.0, 3000.0]
    yRange = [+0.0,+0.2]
    OutFile= OutDir + "Rates.png"
    fig    = mml.plot1D( xRange =xRange , yRange =yRange  , \
                         xTitle =r'$\omega_{ce}t$' , yTitle ="", \
                         FigSize=(5,2.5), FigName=OutFile , \
                         config =config                       )
    fig.addPlot( xAxis=xAxis, yAxis=yAxis1, label="Reconnection Rate" )
    fig.addPlot( xAxis=xAxis, yAxis=yAxis2, label="Merging Rate"      )
    fig.addLegend( loc="lower right" )
    fig.writeFigure()
    
# ------------------------------- #
# --- 実行例 --- #
if ( __name__=='__main__' ):
    import myUtils.myRecvArgs as mra
    args    = mra.myRecvArgs()
    job     = args["job"   ]
    InpDir  = args["jobDir"] + "bin/"
    OutDir  = args["jobDir"] + "png/"
    CnsFile = args["jobDir"] + "dat/" + "constants.dat"
    # exe     = RecRateMaterialDriver( InpDir=InpDir, CnsFile=CnsFile, ksteps=args["ksteps"] )
    exe     = RecRate()
    plt     = RecRatePlotter()
    
    # for kstep in args["ksteps"]:
    #     InpFile = "Field2D{0:08}.bin".format( kstep )
    #     exe     = RecRateMaterial( InpFile=InpFile, InpDir=InpDir, CnsFile=CnsFile )
    #     print( exe )
    # exe = RecRate()
    # RecRatePlotter()
[1;32mkent@gauss[m [1;34m~[m [1;34m$ [m
