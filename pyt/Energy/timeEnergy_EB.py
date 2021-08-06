import numpy               as np
import myStyle.plot1D      as pl1
import myStyle.LoadConfig  as lcf

# --- 運動エネルギーのプロット --- #
def timeEnergy_EB( job   =None, jobDir =None, InpFile=None, OutFile  =None, \
                   OutDir=None, CnsFile=None, config =None, Normalize='Total-Energy' ):
    # -------------------------- #
    # --- [1] 引数チェック   --- #
    # -------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"                .format( config["pic_jobDir"], job )
    if ( InpFile is None ): InpFile = "{0}dat/energy_field.dat".format( jobDir )
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat"   .format( jobDir )
    if ( OutDir  is None ): OutDir  = "{0}png/"                .format( jobDir )
    if ( OutFile is None ): OutFile = OutDir + "timeEnergy_EB.png"
    
    # -------------------------- #
    # --- [2] エネルギー読込 --- #
    # -------------------------- #
    Energy = np.loadtxt( InpFile )
    tAxis  = Energy[:,1]
    if   ( Normalize == "Total-Energy" ):
        WEnergy = np.loadtxt( "{0}dat/energy.dat".format(jobDir) )
        UmaxInv = 1.0 / np.max( WEnergy[:,5] )
    elif ( Normalize == "Total-Field"):
        UmaxInv = 1.0 / np.max(  Energy[:,2] )
    else:
        UmaxInv = 1.0
    Unrm        = Energy[:,3:] * UmaxInv
    Unrm[:,0:4] = Unrm[:,0:4] * 1.e+3
    Uidx        = { "UE1" :0, "UE2":1, "UE3":2, "UET":3, \
                    "UB1" :4, "UB2":5, "UB3":6, "UBT":7  }
    ULabel      = { "UE1":"$U_{e,r}$", "UE2":"$U_{e,t}$", "UE3":"$U_{e,z}$", "UET":"$U_{e,Total}$", \
                    "UB1":"$U_{e,r}$", "UB2":"$U_{e,t}$", "UB3":"$U_{e,z}$", "UBT":"$U_{e,Total}$", \
                    "UTotal":"$U_{Total}$" }
    Utot        = Energy[:,Uidx["UET"]] + Energy[:,Uidx["UBT"]]
    
    # -------------------------- #
    # --- [3]  プロット設定  --- #
    # -------------------------- #
    config["plt_xRange"]      = [0.0, np.max( tAxis ) ]
    config["plt_yRange"]      = [0.0, 1.2             ]
    config["xTitle"]          = "Time"
    config["FigSize"]         = (6,3)
    config["plt_LegLocation"] = "upper right"
    config["plt_LegFontSize"] = 10
    config["plt_LegNColumn"]  = 2

    # -------------------------- #
    # --- [4] プロット       --- #
    # -------------------------- #
    #  -- [4-1] E-Field      --  #
    config["yTitle"]     = "Energy / $10^{-3}$"
    fig    = pl1.plot1D( config=config, FigName=OutFile.replace("EB","E") )
    fig.addPlot( xAxis=tAxis, yAxis=Unrm[:,Uidx["UE1"]], label=ULabel["UE1"] )
    fig.addPlot( xAxis=tAxis, yAxis=Unrm[:,Uidx["UE2"]], label=ULabel["UE2"] )
    fig.addPlot( xAxis=tAxis, yAxis=Unrm[:,Uidx["UE3"]], label=ULabel["UE3"] )
    fig.addPlot( xAxis=tAxis, yAxis=Unrm[:,Uidx["UET"]], label=ULabel["UET"] )
    fig.addLegend()
    fig.writeFigure()
    #  -- [4-2] B-Field      --  #
    config["yTitle"]     = "Energy"
    fig    = pl1.plot1D( config=config, FigName=OutFile.replace("EB","B") )
    fig.addPlot( xAxis=tAxis, yAxis=Unrm[:,Uidx["UB1"]], label=ULabel["UB1"] )
    fig.addPlot( xAxis=tAxis, yAxis=Unrm[:,Uidx["UB2"]], label=ULabel["UB2"] )
    fig.addPlot( xAxis=tAxis, yAxis=Unrm[:,Uidx["UB3"]], label=ULabel["UB3"] )
    fig.addPlot( xAxis=tAxis, yAxis=Unrm[:,Uidx["UBT"]], label=ULabel["UBT"] )
    fig.addLegend()
    fig.writeFigure()

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args      = rar.myRecvArgs()
    OutDir    = "./"
    print( "[timeEnergy_EB] {0} is under processed...".format( args["job"] ), end="" )
    Output    = timeEnergy_EB( job=args["job"], OutDir=OutDir )
    print( "\t [Done]" )
