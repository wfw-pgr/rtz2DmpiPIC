import numpy                   as np
import myStyle.plot1D          as pl1
import myStyle.LoadConfig      as lcf
import mPICsee.TimeUnitConvert as tuc
import myStyle.configSettings  as cfs

def timeEnergy( job   =None, jobDir =None, InpFile=None, OutFile=None, \
                OutDir=None, CnsFile=None, config =None ):
    # -------------------------- #
    # --- [1] 引数チェック   --- #
    # -------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"             .format( config["pic_jobDir"], job )
    if ( InpFile is None ): InpFile = "{0}dat/energy.dat"   .format( jobDir )
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    if ( OutDir  is None ): OutDir  = "{0}png/"             .format( jobDir )
    if ( OutFile is None ): OutFile = OutDir + "timeEnergy.png"
    
    # -------------------------- #
    # --- [2] エネルギー読込 --- #
    # -------------------------- #
    Energy = np.loadtxt( InpFile )
    tAxis  = tuc.TimeUnitConvert( kstep=Energy[:,0], job=job, jobDir=jobDir, display=False, unit="wpi" )
    Emax   = np.max( Energy[:,5] )
    Enrm   = Energy[:,2:6] / Emax
    Eidx   = {"Ek":0, "Ee":1, "Eb":2, "Et":3 }
    
    # -------------------------- #
    # --- [3] プロット設定   --- #
    # -------------------------- #
    cfs.configSettings( config=config, configType="plot1D_lateral" )
    cfs.configSettings( config=config, configType="plot1D_stack" )
    config["plt_xAutoRange"]  = True
    config["plt_yAutoRange"]  = False
    config["plt_xRange"]      = [0.0, 240.]
    config["plt_yRange"]      = [0.0, 1.2 ]
    config["xMajor_Nticks"]   = 7
    config["yMajor_Nticks"]   = 7
    config["xTitle"]          = ""
    config["yTitle"]          = ""
    config["plt_LegLocation"] = "lower right"
    config["plt_LegFontSize"] = 14
    config["plt_LegNColumn"]  = 4
    config["MinimalOut"]      = False
    
    # -------------------------- #
    # --- [4] プロット       --- #
    # -------------------------- #
    fig    = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=tAxis, yAxis=Enrm[:,Eidx["Ek"]], label=r'$U_{Kinetic}$'  )
    fig.addPlot( xAxis=tAxis, yAxis=Enrm[:,Eidx["Ee"]], label=r'$U_{Electric}$' )
    fig.addPlot( xAxis=tAxis, yAxis=Enrm[:,Eidx["Eb"]], label=r'$U_{Magnetic}$' )
    fig.addPlot( xAxis=tAxis, yAxis=Enrm[:,Eidx["Et"]], label=r'$U_{Total}$'    )
    fig.setAxis()
    fig.addLegend()
    fig.writeFigure()

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    OutDir  = None
    print( "[timeEnergy] {0} is under processed...".format( args["job"] ), end="" )
    Output  = timeEnergy( job=args["job"], OutDir=OutDir )
    print( "\t [Done]" )
