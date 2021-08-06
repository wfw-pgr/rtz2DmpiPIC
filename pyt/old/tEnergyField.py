import numpy               as np
import myStyle.myMPL       as mml
import myStyle.LoadConfig  as lcf

# --- 運動エネルギーのプロット --- #
def tEnergyField( job=None, InpFile=None, TotFile=None, OutFile=None, OutDir=None, CnsFile=None ):
    config = lcf.LoadConfig()
    # --- [1] 引数チェック --- #
    if ( InpFile is None ):
        InpFile  = "../../job/{0}/dat/energy_field.dat".format( job )
    if ( TotFile is None ):
        TotFile  = "../../job/{0}/dat/energy.dat".format( job )
    if ( CnsFile is None ):
        CnsFile  = "../../job/{0}/dat/constants.dat".format( job )
    if ( OutFile is None ):
        OutFile  = "tEnergyField.png"
    if ( OutDir  is not None ):
        OutFile  = OutDir + OutFile
    
    # --- [2] エネルギー読込 --- #
    Field  = np.loadtxt( InpFile )
    Total  = np.loadtxt( TotFile )
    UMax   = np.max( Total[:,5]  )
    tAxis  = +   Total[:,1]
    UE1    = + Field[:,2] / UMax * 1e+3
    UE2    = + Field[:,3] / UMax * 1e+3
    UE3    = + Field[:,4] / UMax * 1e+3
    UET    = + Field[:,5] / UMax * 1e+3
    UB1    = + Field[:,6] / UMax
    UB2    = + Field[:,7] / UMax
    UB3    = + Field[:,8] / UMax
    UBT    = + Field[:,9] / UMax
    # --- [3] プロット設定 --- #
    config["plt_xRange"] = [0.0, np.max( tAxis ) ]
    config["xTitle"]     = "Time"
    config["yTitle"]     = "Energy"
    config["FigSize"]    = (6,3)
    config["plt_LegLocation"] = "upper right"
    # --- [4] プロット --- #
    #  -- [4-1] E-Field Out -- #
    config["yTitle"]     = r"Energy / $10^{-3}$"
    config["plt_yRange"] = [0.0, 1.2]
    config["FigName"]    = OutFile.replace( "Field.png", "EField.png" )
    fig                  = mml.myPlot1D( config=config           )
    fig.addPlot( xAxis=tAxis, yAxis=UE1, label=r'$U_{E,r}$'      )
    fig.addPlot( xAxis=tAxis, yAxis=UE2, label=r'$U_{E,t}$'      )
    fig.addPlot( xAxis=tAxis, yAxis=UE3, label=r'$U_{E,z}$'      )
    fig.addPlot( xAxis=tAxis, yAxis=UET, label=r'$U_{E,Total}$'  )
    fig.addLegend()
    fig.writeFigure()
    
    #  -- [4-1] E-Field Out -- #
    config["yTitle"]     = r"Energy"
    config["plt_yRange"] = [0.0, 1.2]
    print( OutFile, config["FigName"] )
    config["FigName"]    = OutFile.replace( "Field.png", "BField.png" )
    print( OutFile, config["FigName"] )
    fig                  = mml.myPlot1D( config=config           )
    fig.addPlot( xAxis=tAxis, yAxis=UB1, label=r'$U_{M,r}$'      )
    fig.addPlot( xAxis=tAxis, yAxis=UB2, label=r'$U_{M,t}$'      )
    fig.addPlot( xAxis=tAxis, yAxis=UB3, label=r'$U_{M,z}$'      )
    fig.addPlot( xAxis=tAxis, yAxis=UBT, label=r'$U_{M,Total}$'  )
    fig.addLegend()
    fig.writeFigure()
    
# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args      = rar.myRecvArgs()
    CnsFile   = args["jobDir"] + "dat/constants.dat"
    OutDir    = "./png/"
    # OutDir    = args["jobDir"] + "png/"
    Output    = tEnergyField( job=args["job"], OutDir=OutDir, CnsFile=CnsFile )
    print( "tEnergy is processed..." )
