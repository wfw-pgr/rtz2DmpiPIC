import numpy                         as np
import myStyle.LoadConfig            as lcf
import mPICsee.FetchPIC              as fpc
import myConvert.arr2dct             as a2d
import myBasicAlgs.cylindricalVolume as clv
import myUtils.progressBar           as pgb

# --- LCFS 内の 密度を全積分して時間発展プロット --- #
def timeDensity( job   =None,  jobDir=None, ksteps =None, OutFile=None, \
                 OutDir=None,  config=None, NewFile=True, \
                 Flag__CallCalculater=True, Flag__CallPlotter=True ):
    # ------------------------------- #
    # --- [1]   引数チェック      --- #
    # ------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( ksteps  is None ): ksteps  = np.arange( config["arg_Iter"] )*config["arg_Step"]+config["arg_Init"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}dat/".format( jobDir )
    if ( OutFile is None ): OutFile = OutDir + "timeDensity.dat"

    # ------------------------------- #
    # --- [2] LCFS 密度合計 算出  --- #
    # ------------------------------- #
    if ( Flag__CallCalculater ):
        Nsteps = len( ksteps )
        if ( NewFile ): ( open( OutFile, "w" ) ).close()    # -- 新規ファイル -- #
        pgb.progressBar( iteration=0, total=Nsteps )        # -- 進捗表示     -- #
        for i,kstep in enumerate( ksteps ):
            ret       = timeDensity_calculater( job=job, jobDir=jobDir, kstep=kstep, config=config )
            pgb.progressBar( iteration=i+1, total=Nsteps )
            with open( OutFile, "ab" ) as f:
                np.savetxt( f, np.array( ret ).reshape( (1,6) ) )
    
    # --------------------------------- #
    # --- [3] 返却 /  描画          --- #
    # --------------------------------- #
    #  -- [3-1] .png 書き出し       --  #
    if ( Flag__CallPlotter ):
        timeDensity_plotter( datFile=OutFile, config=config )


# --- Separatrix 内の 密度を全積分して返却 --- #
def timeDensity_calculater( job=None, jobDir=None, kstep=None, config=None, criterion=0.05 ):
    # ------------------------ #
    # --- [1] 引数チェック --- #
    # ------------------------ #
    if ( config  is None ): config   = lcf.LoadConfig()
    if ( job     is None ): job      = config["pic_job"]
    if ( kstep   is None ): kstep    = config["kstep"  ]
    if ( jobDir  is None ): jobDir   = "{0}{1}/".format( config["pic_jobDir"], job )
    # -------------------------- #
    # --- [2] データ呼び出し --- #
    # -------------------------- #
    Dkeys       = ["psi", "ne0", "ni0"]
    Data        = fpc.FetchPIC( job =job  , jobDir=jobDir, kstep =kstep, \
                                keys=Dkeys, config=config )
    # ------------------------- #
    # --- [3] LCFS 判定     --- #
    # ------------------------- #
    #  -- [3-1] 準備 -- #
    crit_F      = np.max( Data["psi"] ) * criterion
    fMask       = np.ravel( np.where( Data["psi"] > crit_F, 1.0, 0.0 ) )
    volume      = clv.cylindricalVolume( zAxis=Data["xAxis"], rAxis=Data["yAxis"] )
    #  -- [3-2] 計算 -- #
    neTot       = np.dot( np.ravel( Data["ne0"]*volume ), fMask )
    niTot       = np.dot( np.ravel( Data["ni0"]*volume ), fMask )
    volIn       = np.dot( np.ravel( volume )            , fMask )
    neAve       = neTot / volIn
    niAve       = niTot / volIn
    # ------------------------- #
    # --- [4] 返却          --- #
    # ------------------------- #
    return( [ kstep, neTot, niTot, neAve, niAve, volIn ] )


# --- .png 書き出し --- #
def timeDensity_plotter( datFile=None, config=None, OutDir=None ):
    # ------------------------ #
    # --- [1] 引数チェック --- #
    # ------------------------ #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( OutDir  is None ): OutDir  = "./png/"
    if ( datFile is None ): sys.exit( "[timeDensity_plotter] ERROR :: No datFile " )
    Data = np.transpose( np.loadtxt( datFile ) ) 
    
    # ----------------------- #
    # --- [2] プロット準備  -- #
    # ----------------------- #
    #  -- [2-1] ラベリング  -- #
    keys  = ["kstep", "neTot", "niTot", "neAve", "niAve", "volIn" ]
    label = ["kstep", "$\Sigma n_e$", "$\Sigma n_i$", "$n_e^{LCFS}$", "$n_i^{LCFS}$", "Volume" ]
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
    #  -- [3-1] 粒子数 - 積分合計 -- #
    OutFile = OutDir + "nTotal.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["kstep"], yAxis=Data["neTot"], label=label["neTot"]  )
    fig.addPlot( xAxis=Data["kstep"], yAxis=Data["niTot"], label=label["niTot"]  )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    #  -- [3-2] 粒子数 - 平均 -- #
    OutFile = OutDir + "nLCFS.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["kstep"], yAxis=Data["neAve"], label=label["neAve"]  )
    fig.addPlot( xAxis=Data["kstep"], yAxis=Data["niAve"], label=label["niAve"]  )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    #  -- [3-3] LCFS 体積    -- #
    OutFile = OutDir + "volume.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["kstep"], yAxis=Data["volIn"], label=label["volIn"]  )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()


# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    OutDir  = "./"
    print( "[timeDensity] {0} is under processed...".format( args["job"] ), end="\n" )
    Output  = timeDensity( job=args["job"], OutDir=OutDir, ksteps=args["ksteps"] )
    print( "\t [Done]" )
