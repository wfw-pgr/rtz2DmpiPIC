import numpy               as np
import myStyle.plot1D      as pl1
import myStyle.LoadConfig  as lcf

def timeEnergyCompare( job=None, InpFile=None, OutFile=None, OutDir=None, CnsFile=None, config=None ):
    # --- [1] 引数チェック   --- #
    if ( config  is     None ): config = lcf.LoadConfig()
    if ( InpFile is     None ): InpFile  = "../../../job/{0}/dat/energy.dat".format( job )
    if ( CnsFile is     None ): CnsFile  = "../../../job/{0}/dat/constants.dat".format( job )
    if ( OutFile is     None ): OutFile  = "timeEnergyCompare.png"
    if ( OutDir  is not None ): OutFile  = OutDir + OutFile
    # --- [2] エネルギー読込 --- #
    InpFile1     = "../../../job/ctrI02_/dat/energy.dat"
    InpFile2     = "../../../job/ctrO02_/dat/energy.dat"
    Energy1      = np.loadtxt( InpFile1 )
    Energy2      = np.loadtxt( InpFile2 )
    tAxis1       = + Energy1[:,1]
    tAxis2       = + Energy2[:,1]
    Emax         = np.max( Energy1[:,5] )
    idx          = { "Ek":2, "Ee":3, "Eb":4, "Et":5 }
    Enrm1, Enrm2 = Energy1, Energy2
    Enrm1[:,2:6] = Energy1[:,2:6] / Emax
    Enrm2[:,2:6] = Energy2[:,2:6] / Emax
    
    # --- [3] プロット設定 --- #
    config["AutoRange"]       = False
    config["AutoTicks"]       = False
    config["plt_xRange"]      = [0.0, 2000.]
    config["plt_yRange"]      = [0.0, 1.2  ]
    config["axes_x_ticks"]    = [0.,500.,1000.,1500.,2000.]
    config["axes_y_ticks" ]   = [0.,0.25,0.5,0.75,1.]
    config["xTitle"]          = "Time"
    config["yTitle"]          = "Energy"
    config["FigSize"]         = (6,3)
    config["plt_LegLocation"] = "upper right"
    config["plt_LegFontSize"] = 12
    config["plt_LegNColumn"]  = 2
    
    # --- [4] プロット --- #
    fig    = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=tAxis1, yAxis=Enrm1[:,idx["Ek"]], color="crimson"  , label=r'$U_{Kinetic}^{I}$' )
    fig.addPlot( xAxis=tAxis1, yAxis=Enrm1[:,idx["Eb"]], color="magenta"  , label=r'$U_{Magnetic}^{I}$')
    fig.addPlot( xAxis=tAxis2, yAxis=Enrm2[:,idx["Ek"]], color="darkcyan" , label=r'$U_{Kinetic}^{O}$' )
    fig.addPlot( xAxis=tAxis2, yAxis=Enrm2[:,idx["Eb"]], color="royalblue", label=r'$U_{Magnetic}^{O}$')
    fig.addLegend()
    fig.writeFigure()

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args      = rar.myRecvArgs()
    CnsFile   = args["jobDir"] + "dat/constants.dat"
    OutDir    = "./png/"
    Output    = timeEnergyCompare( job=args["job"], OutDir=OutDir, CnsFile=CnsFile )
    print( "timeEnergyCompare is processed..." )
