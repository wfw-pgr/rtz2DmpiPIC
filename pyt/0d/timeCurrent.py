import numpy                  as np
import myStyle.plot1D         as pl1
import myStyle.LoadConfig     as lcf
import myConvert.arr2dct      as a2d
import myStyle.configSettings as cfs

def timeCurrent( job    =None, jobDir=None, InpFile=None, OutFile=None, OutDir=None, \
                 CnsFile=None, config=None ):
    # ----------------- #
    # -- 引数チェック -- #
    # ----------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"             .format( config["pic_jobDir"], job )
    if ( InpFile is None ): InpFile = "{0}dat/Current.dat"  .format( jobDir )
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    if ( OutDir  is None ): OutDir  = "{0}png/"
    if ( OutFile is None ): OutFile = "timeCurrent.png"
    
    # ----------------- #
    # -- データ取得 --  #
    # ----------------- #
    Data  = np.transpose( np.loadtxt( InpFile ) )
    keys = [ "time", \
             "Ie", "Ie95", "IePriv1", "IePriv2", "IeSim", \
             "Ii", "Ii95", "IiPriv1", "IiPriv2", "IiSim", \
             "Ip", "Ip95", "IpPriv1", "IpPriv2", "IpSim"  ]
    Data  = a2d.arr2dct( Data=Data, keys=keys )
    
    # ----------------- #
    # -- コンフィグ  -- #
    # ----------------- #
    cfs.configSettings( config=config, configType="plot1D_def" )
    config["xTitle"]          = "Time"
    config["yTitle"]          = ""
    config["plt_LegFontSize"] = 10
    config["plt_LegNColumn"]  = 1

    # ----------------- #
    # --   プロット  -- #
    # ----------------- #
    # -- Ie95 + Ii95 + Ip95 -- #
    OutFile = OutDir + "Ip95.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["Ie95"] , label=r'$I_{e}^{95}$'  )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["Ii95"] , label=r'$I_{i}^{95}$'  )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["Ip95"] , label=r'$I_{p}^{95}$'  )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()

    # -- Ie   + Ii   + Ip   -- #
    OutFile = OutDir + "Ip.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["Ie"]   , label=r'$I_{e}$'  )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["Ii"]   , label=r'$I_{i}$'  )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["Ip"]   , label=r'$I_{p}$'  )
    fig.addLegend()
    # fig.setAxis()
    fig.writeFigure()

    # -- IeSim, IiSim, IpSim   -- #
    OutFile = OutDir + "Ipsim.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["IeSim"], label=r'$I_{e}$'  )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["IiSim"], label=r'$I_{i}$'  )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["IpSim"], label=r'$I_{p}$'  )
    fig.addLegend()
    # fig.setAxis()
    fig.writeFigure()



if ( __name__=="__main__" ):
    # -- コマンドライン引数 -- #
    import myUtils.myRecvArgs as rar
    args      = rar.myRecvArgs()
    OutDir    = "./png/"
    # -- 実行 -- #
    print( "[timeCurrent] job = {0} is under processing...".format( args["job"] ), end="" )
    timeCurrent( job=args["job"], OutDir=OutDir )
    print( "\t [completed]")
