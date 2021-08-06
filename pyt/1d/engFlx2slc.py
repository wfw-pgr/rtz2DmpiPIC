import numpy                  as np
import myStyle.LoadConfig     as lcf
import mPICsee.FetchPIC       as fpc
import myStyle.plot1D         as pl1
import myStyle.configSettings as cfs
import myConvert.arr2dct      as a2d
import IORoutines.dctarr2dat  as d2d
import IORoutines.LoadConst   as lcn

# ================================================================ #
# ===   engFlx2slc :: エネルギーフラックス  r-1D plot          === #
# ================================================================ #
def engFlx2slc( job   =None, kstep =None, jobDir =None, pngFile=None, pngDir   =None , cnsFile=None, \
                pkeys =None, labels=None, x1Range=None, x2Range=None, writeData=False, \
                config=None  ):
    # ---------------------------------------- #
    # --- [1]    引数チェック              --- #
    # ---------------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( pngDir  is None ): pngDir  = "{0}png/".format( jobDir )
    if ( cnsFile is None ): cnsFile = "{0}dat/constants.dat".format( jobDir )
    if ( pngFile is None ): pngFile = "{0}engFlx_{1}__{2:08}.png".format( pngDir, job, kstep )
    if ( pkeys   is None ): pkeys   = ["Fix","ni","uix","SFix"]
    if ( labels  is None ): labels  = [None]*len( pkeys )
    Flag_qthermal = True
    pltAxis = "x2"
    Fkeys   = ["ne","uex","uey","uez","pexx","peyy","pezz","pexy", "peyz", "pexz", \
               "ni","uix","uiy","uiz","pixx","piyy","pizz","pixy", "piyz", "pixz", ]
    if ( Flag_qthermal ):
        Fkeys = Fkeys + [ "qvex", "qvix" ]
        config["pic_nMoments"] = 16
    x1Range = [-1.0,1.0]
    x2Range = None
    # ---------------------------------------- #
    # --- [2] データ呼び出し               --- #
    # ---------------------------------------- #
    #  -- 軸の選択 / 平均化  --  #
    if ( pltAxis=="x2" ): # x2 :: LJ, r
        meanAxis = 1
        if ( x1Range is None ): x1Range = config["Data_avgRange"]
        if ( x2Range is None ): x2Range = config["Data_x1Range" ]
    else:
        import sys
        sys.exit( "[pic2slc]  ERROR!! :: plt_Axis??" )
    #  -- データ読み出し     --  #
    config["cmp_LinearFilt"]  = 0.05
    config["cmp_nFilter"]     = 100
    const        = lcn.LoadConst( FileName=cnsFile )
    Data         = fpc.FetchPIC( job =job  , jobDir=jobDir, kstep  =kstep, \
                                 keys=Fkeys, config=config, x1Range=x1Range, x2Range=x2Range )
    xAxis        = Data["xAxis"] if ( pltAxis == "x1" ) else Data["yAxis"]
    for i, key in enumerate(Fkeys): Data[key] = np.ravel( np.mean( Data[key], axis=meanAxis ) )
    print( const["rmi"], const["rme"] )
    Flx          = {}
    Flx["FEke"]  = 0.5 * const["rme"] * Data["ne"] * ( Data["uex"]**2 + Data["uey"]**2 + Data["uez"]**2 ) * Data["uex"]
    Flx["FEki"]  = 0.5 * const["rmi"] * Data["ni"] * ( Data["uix"]**2 + Data["uiy"]**2 + Data["uiz"]**2 ) * Data["uix"]
    Flx["Fpue"]  = 1.5 * ( Data["pexx"] + Data["peyy"] + Data["pezz"] ) / 3.0 * Data["uex"] + Data["pexx"]*Data["uex"]
    Flx["Fpui"]  = 1.5 * ( Data["pixx"] + Data["piyy"] + Data["pizz"] ) / 3.0 * Data["uix"] + Data["pixx"]*Data["uix"]
    Flx["Fnde"]  = ( Data["pexy"]*Data["uey"] + Data["pexz"]*Data["uez"] )
    Flx["Fndi"]  = ( Data["pixy"]*Data["uiy"] + Data["pixz"]*Data["uiz"] ) 
    if   ( Flag_qthermal ):
        Flx["Fqte"]  = Data["qvex"]
        Flx["Fqti"]  = Data["qvix"]
    else:
        Flx["Fqte"]  = 0.0 * xAxis
        Flx["Fqti"]  = 0.0 * xAxis
    ikeys        = ["FEki","Fpui","Fndi","Fqti"]
    ekeys        = ["FEke","Fpue","Fnde","Fqte"]
    ilabels      = a2d.arr2dct( Data=ikeys, keys=ikeys )
    elabels      = a2d.arr2dct( Data=ekeys, keys=ekeys )
    FluxDensity  = True
    if ( not( FluxDensity ) ):
        sfcArea      = 2.0 * np.pi * xAxis
        print( sfcArea.shape )
        for key in ikeys: Flx[key] = Flx[key] * sfcArea
        for key in ekeys: Flx[key] = Flx[key] * sfcArea

    # ---------------------------------------- #
    # --- [3] コンフィグ設定               --- #
    # ---------------------------------------- #
    cfs.configSettings( configType="plot1D_def"    , config=config )
    cfs.configSettings( configType="plot1D_lateral", config=config )
    cfs.configSettings( configType="plot1D_stack"  , config=config )
    cfs.configSettings( configType="NoAxis"        , config=config )
    config["FigSize"]        = (6,6)
    config["plt_yAutoRange"] = False
    # ---------------------------------------- #
    # --- [4] 流束密度/流速/密度           --- #
    # ---------------------------------------- #
    pngFile_ = pngFile.replace( "__", "engFlx_i" )
    # config["plt_yRange"]     = [-2.0,2.0]
    config["plt_yRange"]     = [-0.05,0.05]
    fig      = pl1.plot1D( FigName=pngFile_, config=config )
    for key in ikeys: fig.addPlot( xAxis=xAxis, yAxis=Flx[key], label=ilabels[key] )
    fig.setAxis()
    fig.writeFigure()
    # ---------------------------------------- #
    # --- [5] 流束 ( 表面積考慮 )          --- #
    # ---------------------------------------- #
    pngFile_ = pngFile.replace( "__", "engFlx_e" )
    # config["plt_yRange"]     = [-2.0,2.0]
    config["plt_yRange"]     = [-0.05,0.05]
    fig          = pl1.plot1D( FigName=pngFile_, config=config )
    for key in ekeys: fig.addPlot( xAxis=xAxis, yAxis=Flx[key], label=elabels[key] )
    fig.setAxis()
    fig.writeFigure()
    # ---------------------------------------- #
    # --- [6] データ書出                   --- #
    # ---------------------------------------- #
    if ( writeData ):
        datFile = pngFile.replace( ".png", ".dat" )
        d2d.dctarr2dat( Data=Flx, keys=ekeys+ikeys, datFile=datFile )

# ================================================================ #
# ===  実行部                                                  === #
# ================================================================ #
if ( __name__=="__main__" ):
    import myUtils.genArgs as gar
    args    = gar.genArgs()
    for kstep in args["ksteps"]:
        engFlx2slc( job    =args["job"]   , jobDir=args["jobDir"], kstep =kstep, \
                    pngDir =args["pngDir"], pkeys =args["key"   ] )
