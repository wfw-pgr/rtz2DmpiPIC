import numpy                  as np
import myStyle.plot1D         as pl1
import myStyle.LoadConfig     as lcf
import myConvert.arr2dct      as a2d
import myStyle.configSettings as cfs

# ---- Momentum トロイダル角運動量 を時間発展でプロット ---- #
def timeMomentum( job    =None, jobDir=None, InpFile=None, OutFile=None, OutDir=None, \
                  CnsFile=None, config=None ):
    # ----------------- #
    # -- 引数チェック -- #
    # ----------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"                   .format( config["pic_jobDir"], job )
    if ( InpFile is None ): InpFile = "{0}dat/AngularMomentum.dat".format( jobDir )
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat"      .format( jobDir )
    if ( OutDir  is None ): OutDir  = "{0}png/"
    if ( OutFile is None ): OutFile = "timeCurrent.png"
    
    # ---------------- #
    # -- データ取得 -- #
    # ---------------- #
    Data   = np.transpose( np.loadtxt( InpFile ) )
    keys   = [ "time", \
               "aMe", "aMe95", "aMePriv1", "aMePriv2", "aMeSim", \
               "aMi", "aMi95", "aMiPriv1", "aMiPriv2", "aMiSim"  ]
    Data   = a2d.arr2dct( Data=Data, keys=keys )
    MLabel = { "time"    :"Time", \
               "aMe"     :"$m_{e}rv_{e,\phi}$"        , "aMi"     :"$m_{i}rv_{i,\phi}$"        , \
               "aMe95"   :"$m_{e}rv_{e,\phi}^{LCFS}$" , "aMi95"   :"$m_{i}rv_{i,\phi}^{LCFS}$" , \
               "aMePriv1":"$m_{e}rv_{e,\phi}^{Sph.1}$", "aMiPriv1":"$m_{i}rv_{i,\phi}^{Sph.1}$", \
               "aMePriv2":"$m_{e}rv_{e,\phi}^{Sph.2}$", "aMiPriv2":"$m_{i}rv_{i,\phi}^{Sph.2}$", \
               "aMeSim"  :"$m_{e}rv_{e,\phi}^{Total}$", "aMiSim"  :"$m_{i}rv_{i,\phi}^{Total}$"  }
    
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
    # -- LCFS -- #
    OutFile = OutDir + "AnglarMomentum.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["aMe"], label=MLabel["aMe"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["aMi"], label=MLabel["aMi"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    # -- LCFS 95 -- #
    OutFile = OutDir + "AnglarMomentum_95.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["aMe95"], label=MLabel["aMe95"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["aMi95"], label=MLabel["aMi95"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    # -- LCFS Private -- #
    OutFile = OutDir + "AnglarMomentum_Priv.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["aMePriv1"], label=MLabel["aMePriv1"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["aMiPriv1"], label=MLabel["aMiPriv1"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["aMePriv2"], label=MLabel["aMePriv2"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["aMiPriv2"], label=MLabel["aMiPriv2"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    # -- LCFS Entire Sim. -- #
    OutFile = OutDir + "AnglarMomentum_Sim.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["aMeSim"], label=MLabel["aMeSim"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["aMiSim"], label=MLabel["aMiSim"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()



# ----------------------------------------------------- #
if ( __name__=="__main__" ):
    # -- コマンドライン引数 -- #
    import myUtils.myRecvArgs as rar
    args      = rar.myRecvArgs()
    OutDir    = "./png/"
    # -- 実行 -- #
    print( "[timeMomentum] {0} is under processing...".format( args["job"] ), end="" )
    timeMomentum( job=args["job"], OutDir=OutDir )
    print( "\t [completed]")
