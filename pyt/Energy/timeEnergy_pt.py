import numpy               as np
import myStyle.plot1D      as pl1
import myStyle.LoadConfig  as lcf

# --- 運動エネルギーのプロット --- #
def timeEnergy_pt( job   =None, jobDir =None, InpFile=None, OutFile  =None, \
                   OutDir=None, CnsFile=None, config =None, Normalize='Total-Energy' ):
    # ----------------------- #
    # --- [1] 引数チェック --- #
    # ----------------------- #
    if ( config  is None ): config   = lcf.LoadConfig()
    if ( job     is None ): job      = config["pic_job"]
    if ( jobDir  is None ): jobDir   = "{0}{1}/"                  .format( config["pic_jobDir"], job )
    if ( InpFile is None ): InpFile  = "{0}dat/energy_kinetic.dat".format( jobDir )
    if ( CnsFile is None ): CnsFile  = "{0}dat/constants.dat"     .format( jobDir )
    if ( OutDir  is None ): OutDir   = "{0}png/"                  .format( jobDir )
    if ( OutFile is None ): OutFile  = OutDir + "timeEnergy_pt.png"
    
    # ------------------------- #
    # --- [2] エネルギー読込 --- #
    # ------------------------- #
    Energy = np.loadtxt( InpFile )
    tAxis  = Energy[:,1]
    if   ( Normalize == "Total-Energy" ):
        WEnergy = np.loadtxt( "{0}dat/energy.dat".format(jobDir) )
        Kmax    = np.max( WEnergy[:,5] )
    elif ( Normalize == "Total-Kinetic"):
        Kmax    = np.max( Energy[:,2] )
    else:
        Kmax    = 1.0
    Knrm   = Energy[:,2:] / Kmax
    Kidx   = { "Kex" :1, "Key":2, "Kez":3, "Ke":4, \
               "Kix" :5, "Kiy":6, "Kiz":7, "Ki":8, "Ktotal":0 }
    KLabel = { "Kex":"$K_{e,r}$", "Key":"$K_{e,t}$", "Kez":"$K_{e,z}$", "Ke":"$K_{e,Total}$", \
               "Kix":"$K_{i,r}$", "Kiy":"$K_{i,t}$", "Kiz":"$K_{i,z}$", "Ki":"$K_{i,Total}$", "Ktotal":"$K_{e+i}$"}
    
    # ------------------------- #
    # --- [3]  プロット設定  --- #
    # ------------------------- #
    config["plt_xRange"] = [0.0, np.max( tAxis ) ]
    config["plt_yRange"] = [0.0, 1.2             ]
    config["xTitle"]     = "Time"
    config["yTitle"]     = "Energy"
    config["FigSize"]    = (6,3)
    config["plt_LegLocation"] = "upper left"
    config["plt_LegFontSize"] = 12
    config["plt_LegNColumn"]  = 4

    # --- [4] プロット --- #
    fig    = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=tAxis, yAxis=Knrm[:,Kidx["Kex"]]   , label=KLabel["Kex"]     )
    fig.addPlot( xAxis=tAxis, yAxis=Knrm[:,Kidx["Key"]]   , label=KLabel["Key"]     )
    fig.addPlot( xAxis=tAxis, yAxis=Knrm[:,Kidx["Kez"]]   , label=KLabel["Kez"]     )
    fig.addPlot( xAxis=tAxis, yAxis=Knrm[:,Kidx["Ke" ]]   , label=KLabel["Ke" ]     )
    fig.addPlot( xAxis=tAxis, yAxis=Knrm[:,Kidx["Kix"]]   , label=KLabel["Kix"]     )
    fig.addPlot( xAxis=tAxis, yAxis=Knrm[:,Kidx["Kiy"]]   , label=KLabel["Kiy"]     )
    fig.addPlot( xAxis=tAxis, yAxis=Knrm[:,Kidx["Kiz"]]   , label=KLabel["Kiz"]     )
    fig.addPlot( xAxis=tAxis, yAxis=Knrm[:,Kidx["Ki" ]]   , label=KLabel["Ki" ]     )
    fig.addPlot( xAxis=tAxis, yAxis=Knrm[:,Kidx["Ktotal"]], label=KLabel["Ktotal"]  )
    fig.addLegend()
    fig.writeFigure()

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args      = rar.myRecvArgs()
    OutDir    = "./"
    print( "[timeEnergy_pt] {0} is under processed...".format( args["job"] ), end="" )
    Output    = timeEnergy_pt( job=args["job"], OutDir=OutDir )
    print( "\t [Done]" )
