import numpy                  as np
import myStyle.plot1D         as pl1
import myStyle.LoadConfig     as lcf
import myConvert.arr2dct      as a2d
import myStyle.configSettings as cfs

# ---- Flux ( phi, psi ) を時間発展でプロット ---- #
def timeFlux( job    =None, jobDir=None, InpFile=None, OutFile=None, OutDir=None, \
              CnsFile=None, config=None ):
    # ----------------- #
    # -- 引数チェック -- #
    # ----------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"                .format( config["pic_jobDir"], job )
    if ( InpFile is None ): InpFile = "{0}dat/MagneticFlux.dat".format( jobDir )
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat"   .format( jobDir )
    if ( OutDir  is None ): OutDir  = "{0}png/"                .format( jobDir )
    if ( OutFile is None ): OutFile = "timeFlux.png"
    
    # ---------------- #
    # -- データ取得 -- #
    # ---------------- #
    Data   = np.transpose( np.loadtxt( InpFile ) )
    keys   = [ "time" , "psiO1", "psiO2" , "psiFX"   , "psiNX"   , "psiMin", \
               "rOpt1", "zOpt1", "rOpt2" , "zOpt2"   , "rXpt"    , "zXpt"  , \
               "phi95", "phi"  , "phiSim", "phiPriv1", "phiPriv2", "mergR" ]
    Data   = a2d.arr2dct( Data=Data, keys=keys )
    PLabel = { "time"    :"$Time$"         , "psiO1"   :"$\psi_{O}^{Sph.1}$", "psiO2" :"$\psi_{O}^{Sph.1}$", \
               "psiFX"   :"$\psi_{X}$"     , "psiNX"   :"$\psi_{X}^{norm}$" , "psiMin":"$\psi_{Min}$"      , \
               "rOpt1"   :"$R_{O}^{Sph.1}$", "rOpt2"   :"$R_{O}^{Sph.2}$"   , "rXpt"  :"$R_{X}$"           , \
               "zOpt1"   :"$Z_{O}^{Sph.1}$", "zOpt2"   :"$Z_{O}^{Sph.2}$"   , "zXpt"  :"$Z_{X}$"           , \
               "phi95"   :"$\phi^{LCFS95}$", "phi"     :"$\phi^{LCFS}$"     , "phiSim":"$\phi^{Sim.}$"     , \
               "phiPriv1":"$\phi^{Sph.1}$" , "phiPriv2":"$\phi^{Sph.1}$"    , "mergR" :"$\eta_{Merg}$"     }
    
    # ---------------- #
    # -- コンフィグ -- #
    # ---------------- #
    cfs.configSettings( config=config, configType="plot1D_def" )
    config["xTitle"]          = "Time"
    config["yTitle"]          = "Energy"
    config["plt_LegFontSize"] = 10
    config["plt_LegNColumn"]  = 1

    # ---------------- #
    # -- プロット   -- #
    # ---------------- #
    # -- psi O1,O2,X -- #
    OutFile = OutDir + "MagneticFlux.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["psiO1"], label=PLabel["psiO1"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["psiO2"], label=PLabel["psiO2"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["psiFX"], label=PLabel["psiFX"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    # -- Ro, Zo -- #
    OutFile = OutDir + "Opoint.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["rOpt1"], label=PLabel["rOpt1"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["zOpt1"], label=PLabel["zOpt1"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["rOpt2"], label=PLabel["rOpt2"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["zOpt2"], label=PLabel["zOpt2"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    # -- Rx, Zx -- #
    OutFile = OutDir + "Xpoint.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["rXpt" ], label=PLabel["rXpt" ] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["zXpt" ], label=PLabel["zXpt" ] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    # -- Merging-Rate -- #
    OutFile = OutDir + "MergingRate.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["mergR" ], label=PLabel["mergR" ] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    # -- phi -- #
    OutFile = OutDir + "ToroidalFlux.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["phi95"   ], label=PLabel["phi95"   ] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["phi"     ], label=PLabel["phi"     ] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["phiSim"  ], label=PLabel["phiSim"  ] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["phiPriv1"], label=PLabel["phiPriv1"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["phiPriv2"], label=PLabel["phiPriv2"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()



# ----------------------------------------------------- #
if ( __name__=="__main__" ):
    # -- コマンドライン引数 -- #
    import myUtils.myRecvArgs as rar
    args      = rar.myRecvArgs()
    OutDir    = None
    # -- 実行 -- #
    print( "[timeFlux] {0} is under processing...".format( args["job"] ), end="" )
    timeFlux( job=args["job"], OutDir=OutDir )
    print( "\t [completed]")
