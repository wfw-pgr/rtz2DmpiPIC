import numpy                  as np
import myStyle.plot1D         as pl1
import myStyle.LoadConfig     as lcf
import myUtils.LoadConst      as lcn
import myConvert.arr2dct      as a2d
import myStyle.configSettings as cfs

# --- エネルギーを0次元解析 ( 時間-エネルギー ) --- #
def timeEnergy( job=None, jobDir=None, InpFile=None, OutFile=None, OutDir=None, CnsFile=None, config=None ):
    # ---------------------------------------- #
    # --- [1] 引数チェック                 --- #
    # ---------------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"             .format( config["pic_jobDir"], job )
    if ( InpFile is None ): InpFile = "{0}dat/energy.dat"   .format( jobDir )
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    if ( OutDir  is None ): OutDir  = "{0}png/"             .format( jobDir )
    if ( OutFile is None ): OutFile = OutDir + "timeEnergy.png"
    
    # ---------------------------------------- #
    # --- [2] エネルギー読込               --- #
    # ---------------------------------------- #
    Const  = lcn.LoadConst( InpFile=CnsFile )
    Data   = np.transpose( np.loadtxt( InpFile ) )
    tAxis  = Data[1,:] / Const["mr"]
    Unrm   = Data[2:6,:] / ( np.max( Data[5,:] ) )
    keys   = ["Uk", "Ue", "Ub", "Ut" ]
    Data   = a2d.arr2dct( Data=Unrm, keys=keys )

    # ---------------------------------------- #
    # --- [3] プロット設定                 --- #
    # ---------------------------------------- #
    cfs.configSettings( config=config, configType="plot1D_def" )
    config["plt_yAutoRange"]  = False
    config["plt_yRange"]      = [0.0,1.2]
    config["xTitle"]          = "Time"
    config["yTitle"]          = "Energy"
    config["plt_LegNColumn"]  = 2
    config["plt_LegFontSize"] = 16

    # ---------------------------------------- #
    # --- [4] プロット設定                 --- #
    # ---------------------------------------- #
    fig    = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=tAxis, yAxis=Data["Uk"], label=r'$U_{Kinetic}$'  )
    fig.addPlot( xAxis=tAxis, yAxis=Data["Ue"], label=r'$U_{Electric}$' )
    fig.addPlot( xAxis=tAxis, yAxis=Data["Ub"], label=r'$U_{Magnetic}$' )
    fig.addPlot( xAxis=tAxis, yAxis=Data["Ut"], label=r'$U_{Total}$'    )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()

# --------------------------------------- #
if ( __name__=="__main__" ):
    # -- コマンドライン引数 -- #
    import myUtils.myRecvArgs as rar
    args      = rar.myRecvArgs()
    OutDir    = None
    # -- 実行 -- #
    timeEnergy( job=args["job"], OutDir=OutDir )
