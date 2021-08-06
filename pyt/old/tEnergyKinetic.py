import numpy               as np
import myStyle.myMPL       as mml
import myStyle.LoadConfig  as lcf

# --- 運動エネルギーのプロット --- #
def tEnergyKinetic( job=None, InpFile=None, TotFile=None, OutFile=None, OutDir=None, CnsFile=None ):
    config = lcf.LoadConfig()
    # --- [1] 引数チェック --- #
    if ( InpFile is None ):
        InpFile  = "../../job/{0}/dat/energy_kinetic.dat".format( job )
    if ( TotFile is None ):
        TotFile  = "../../job/{0}/dat/energy.dat".format( job )
    if ( CnsFile is None ):
        CnsFile  = "../../job/{0}/dat/constants.dat".format( job )
    if ( OutFile is None ):
        OutFile  = "tEnergyKinetic.png"
    if ( OutDir  is not None ):
        OutFile  = OutDir + OutFile
    
    # --- [2] エネルギー読込 --- #
    kinetic= np.loadtxt( InpFile )
    Total  = np.loadtxt( TotFile )
    UMax   = np.max( Total[:,5]  )
    tAxis  = +   Total[:,1]
    UkT    = + kinetic[:,1] / UMax
    Ue1    = + kinetic[:,2] / UMax
    Ue2    = + kinetic[:,3] / UMax
    Ue3    = + kinetic[:,4] / UMax
    UeT    = + kinetic[:,5] / UMax
    Ui1    = + kinetic[:,6] / UMax
    Ui2    = + kinetic[:,7] / UMax
    Ui3    = + kinetic[:,8] / UMax
    UiT    = + kinetic[:,9] / UMax
    # --- [3] プロット設定 --- #
    config["plt_xRange"] = [0.0, np.max( tAxis ) ]
    config["plt_yRange"] = [0.0, 1.2             ]
    config["xTitle"]     = "Time"
    config["yTitle"]     = "Energy"
    config["FigSize"]    = (6,3)
    config["plt_LegLocation"] = "upper right"
    # --- [4] プロット --- #
    fig    = mml.myPlot1D( config=config, FigName=OutFile        )
    fig.addPlot( xAxis=tAxis, yAxis=UkT, label=r'$U_{Kinetic}$'  )
    fig.addPlot( xAxis=tAxis, yAxis=Ue1, label=r'$U_{e,r}$'      )
    fig.addPlot( xAxis=tAxis, yAxis=Ue2, label=r'$U_{e,t}$'      )
    fig.addPlot( xAxis=tAxis, yAxis=Ue3, label=r'$U_{e,z}$'      )
    fig.addPlot( xAxis=tAxis, yAxis=UeT, label=r'$U_{e,Total}$'  )
    fig.addPlot( xAxis=tAxis, yAxis=Ui1, label=r'$U_{i,r}$'      )
    fig.addPlot( xAxis=tAxis, yAxis=Ui2, label=r'$U_{i,t}$'      )
    fig.addPlot( xAxis=tAxis, yAxis=Ui3, label=r'$U_{i,z}$'      )
    fig.addPlot( xAxis=tAxis, yAxis=UiT, label=r'$U_{i,Total}$'  )
    fig.addLegend()
    fig.writeFigure()

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args      = rar.myRecvArgs()
    CnsFile   = args["jobDir"] + "dat/constants.dat"
    OutDir    = "./"
    # OutDir    = args["jobDir"] + "png/"
    Output    = tEnergyKinetic( job=args["job"], OutDir=OutDir, CnsFile=CnsFile )
    print( "tEnergy is processed..." )
