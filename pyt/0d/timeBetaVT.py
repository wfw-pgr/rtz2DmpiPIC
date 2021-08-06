import sys
import numpy                         as np
import myStyle.LoadConfig            as lcf
import mPICsee.FetchPIC              as fpc
import myConvert.arr2dct             as a2d
import myBasicAlgs.cylindricalVolume as clv
import myUtils.progressBar           as pgb
import myUtils.LoadConst             as lcn

# --- LCFS 内の 密度を全積分して時間発展プロット --- #
def timeBetaVT( job   =None,  jobDir=None, ksteps =None, OutFile=None, \
                OutDir=None,  config=None, NewFile=True, CnsFile=None, \
                Flag__CallCalculater=True, Flag__CallPlotter=True ):
    # ------------------------------- #
    # --- [1]   引数チェック      --- #
    # ------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( ksteps  is None ): ksteps  = np.arange( config["arg_Iter"] )*config["arg_Step"]+config["arg_Init"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"             .format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}dat/"             .format( jobDir )
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    if ( OutFile is None ): OutFile = OutDir + "timeBetaVT.dat"

    # ------------------------------- #
    # --- [2] LCFS 密度合計 算出  --- #
    # ------------------------------- #
    if ( Flag__CallCalculater ):
        Nsteps = len( ksteps )
        if ( NewFile ): ( open( OutFile, "w" ) ).close()    # -- 新規ファイル -- #
        pgb.progressBar( iteration=0, total=Nsteps )        # -- 進捗表示     -- #
        for i,kstep in enumerate( ksteps ):
            ret       = timeBetaVT_calculater( job   =job   , jobDir =jobDir, kstep=kstep, \
                                               config=config, CnsFile=CnsFile )
            pgb.progressBar( iteration=i+1, total=Nsteps )
            with open( OutFile, "ab" ) as f:
                np.savetxt( f, np.array( ret ).reshape( (1,14) ) )
    
    # --------------------------------- #
    # --- [3] 返却 /  描画          --- #
    # --------------------------------- #
    #  -- [3-1] .png 書き出し       --  #
    if ( Flag__CallPlotter ):
        timeBetaVT_plotter( datFile=OutFile, config=config )


# --- Separatrix 内の 密度を全積分して返却 --- #
def timeBetaVT_calculater( job   =None, jobDir   =None, kstep=None, CnsFile=None, \
                           config=None, criterion=0.05, gamma=5.0/3.0 ):
    # ------------------------ #
    # --- [1] 引数チェック --- #
    # ------------------------ #
    if ( config  is None ): config   = lcf.LoadConfig()
    if ( job     is None ): job      = config["pic_job"]
    if ( kstep   is None ): kstep    = config["kstep"  ]
    if ( jobDir  is None ): jobDir   = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( CnsFile is None ): CnsFile  = "{0}dat/constants.dat"
    # -------------------------- #
    # --- [2] データ呼び出し --- #
    # -------------------------- #
    const       = lcn.LoadConst( InpFile=CnsFile )
    Dkeys       = ["psi", "Bx", "By", "Bz", "pexx", "peyy", "pezz", "pixx", "piyy", "pizz"]
    Data        = fpc.FetchPIC( job =job  , jobDir=jobDir, kstep =kstep, \
                                keys=Dkeys, config=config )
    time        = const["dt"] * float( kstep )
    
    # ------------------------- #
    # --- [3] LCFS 判定     --- #
    # ------------------------- #
    #  -- [3-1] 準備 -- #
    crit_F      = np.max( Data["psi"] ) * criterion
    fMask       = np.ravel( np.where( Data["psi"] > crit_F, 1.0, 0.0 ) )
    volume      = clv.cylindricalVolume( zAxis=Data["xAxis"], rAxis=Data["yAxis"] )
    vol         = np.ravel( volume )
    #  -- [3-2] 計算 -- #
    coef        = 0.5 / const["wpewce"]**2   # -- 1/2 valfe**2 -- #
    Pmag        = np.ravel( coef*( Data["Bx"]**2 + Data["By"]**2 + Data["Bz"]**2 )*volume )
    PBtf        = np.ravel( coef*(                 Data["By"]**2                 )*volume )
    Pe          = np.ravel(      ( Data["pexx"]  + Data["peyy"]  + Data["pezz"]  )*volume )
    Pi          = np.ravel(      ( Data["pixx"]  + Data["piyy"]  + Data["pizz"]  )*volume )
    PmTot       = np.dot( Pmag, fMask )
    BtTot       = np.dot( PBtf, fMask )
    PeTot       = np.dot( Pe  , fMask )
    PiTot       = np.dot( Pi  , fMask )
    PtTot       = PeTot + PiTot
    volIn       = np.dot( vol , fMask )
    PmAve       = PmTot / volIn
    BtAve       = BtTot / volIn
    PeAve       = PeTot / volIn
    PiAve       = PiTot / volIn
    PtAve       = PtTot / volIn
    BetaV       = PtAve / PmAve
    BetaT       = PtAve / BtAve
    # ------------------------- #
    # --- [4] 返却          --- #
    # ------------------------- #
    return( [ kstep, time , BetaV, BetaT, \
              PmTot, BtTot, PeTot, PiTot, PtTot, \
              PmAve, BtAve, PeAve, PiAve, PtAve  ] )


# --- .png 書き出し --- #
def timeBetaVT_plotter( datFile=None, config=None, OutDir=None ):
    # ------------------------ #
    # --- [1] 引数チェック --- #
    # ------------------------ #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( OutDir  is None ): OutDir  = "./png/"
    if ( datFile is None ): sys.exit( "[timeBetaVT_plotter] ERROR :: No datFile " )
    Data = np.transpose( np.loadtxt( datFile ) ) 
    
    # ------------------------ #
    # --- [2] プロット準備  -- #
    # ------------------------ #
    #  -- [2-1] ラベリング  -- #
    keys  = ["kstep", "time" , "BetaV", "BetaT", \
             "PmTot", "BtTot", "PeTot", "PiTot", "PtTot", \
             "PmAve", "BtAve", "PeAve", "PiAve", "PtAve"  ]
    label = ["kstep", "Time" , "$Beta_{vol}$" , "$Beta_{tor}$", \
             "$\Sigma B^2/2$", "$\Sigma B_t^2/2$", "$\Sigma p_e$", "$\Sigma p_i$", "$\Sigma p$", \
             "$<B^2/2>$", "$<B_t^2/2>$", "$<p_e>$", "$<p_i>$", "$<p>$"]
    Data  = a2d.arr2dct( Data=Data , keys=keys )
    label = a2d.arr2dct( Data=label, keys=keys )
    #  -- [2-2] コンフィグ  -- #
    import myStyle.plot1D         as pl1
    import myStyle.configSettings as cfs
    cfs.configSettings( config=config, configType="plot1D_def" )
    config["xTitle"]          = "Time"
    config["yTitle"]          = ""
    config["plt_LegFontSize"] = 10
    config["plt_LegNColumn"]  = 1

    # ----------------------- #
    # --- [3] プロット    --- #
    # ----------------------- #
    #  -- [3-1] 圧力 -- 単純合計 -- #
    OutFile = OutDir + "sumPBpres.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["PeTot"], label=label["PeTot"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["PiTot"], label=label["PiTot"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["PtTot"], label=label["PtTot"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["PmTot"], label=label["PmTot"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["BtTot"], label=label["BtTot"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    #  -- [3-2] 圧力 --   平均   -- #
    OutFile = OutDir + "avePBpres.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["PeAve"], label=label["PeAve"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["PiAve"], label=label["PiAve"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["PtAve"], label=label["PtAve"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["PmAve"], label=label["PmAve"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["BtAve"], label=label["BtAve"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    #  -- [3-3] 体積平均ベータ値 -- #
    OutFile = OutDir + "BetaTor.png"
    fig     = pl1.plot1D( xAxis=Data["time"], yAxis=Data["BetaT"], label=label["BetaT"], config=config, FigName=OutFile )
    OutFile = OutDir + "BetaVol.png"
    fig     = pl1.plot1D( xAxis=Data["time"], yAxis=Data["BetaV"], label=label["BetaV"], config=config, FigName=OutFile )


# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    OutDir  = "./dat/"
    print( "[timeBetaVT] {0} is under processed...".format( args["job"] ), end="\n" )
    Output  = timeBetaVT( job=args["job"], OutDir=OutDir, ksteps=args["ksteps"], \
                          Flag__CallCalculater=args["calcMode"], Flag__CallPlotter=args["plotMode"] )
    print( "\t [Done]" )
