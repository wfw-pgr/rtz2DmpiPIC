import numpy                   as np
import myStyle.LoadConfig      as lcf
import mPICsee.TimeUnitConvert as tuc
import myStyle.plot1D          as pl1
import myStyle.configSettings  as cfs
import myConvert.arr2dct       as a2d

# ======================================== #
# ===  X-point のr位置 ( rXpt ) を図示 === #
# ======================================== #
def rXptPlotter( jobBase=None, pngDir=None, pngFile=None, config=None ):
    # ---------------------------------------- #
    # --- [1]   引数チェック               --- #
    # ---------------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( jobBase is None ): jobBase = config["pic_jobDir"]
    if ( pngDir  is None ): pngDir  = "./"
    if ( pngFile is None ): pngFile = "{0}rXptPlot.png".format( pngDir )
    # ------------------ #
    # -- 対象ジョブ名 -- #
    # ------------------ #
    jobs  = ["CtrO21_01", "CtrI21_01", "Co-21_01"   ]
    Legs  = ["Case-O"   , "Case-I"   , "Co-Helicity"]
    Legs  = a2d.arr2dct( Data=Legs, keys=jobs )

    # ---------------------------------------- #
    # --- [2] データ読込                   --- #
    # ---------------------------------------- #
    Data  = []
    for job in jobs:
        jobDir    = "{0}{1}/".format( jobBase, job )
        datFile   = "{0}dat/RateAnalysis.dat".format( jobDir )
        hData     = np.transpose( np.loadtxt( datFile, comments="#" ) )
        with open( datFile, "r" ) as f:
            items = ( ( f.readline() ).replace("#","") ).split()
        hData     = a2d.arr2dct( Data=hData, keys=items )
        hData["job"]   = job
        hData["ptime"] = tuc.TimeUnitConvert( kstep=hData["fstep"], unit="wpi", job=job, jobDir=jobDir )
        Data.append( hData )
        
    # ---------------------------------------- #
    # --- [3] コンフィグ設定               --- #
    # ---------------------------------------- #
    cfs.configSettings( config=config, configType="plot1D_def"  )
    cfs.configSettings( config=config, configType="FilterClear" )
    # cfs.configSettings( config=config, configType="NoAxis"      )
    config["FigSize"]        = (8,3)
    config["xTitle"]         = ""
    config["yTitle"]         = ""
    config["plt_xAutoRange"] = False
    config["plt_xRange"]     = [0.0,120.0]
    config["plt_yAutoRange"] = False
    config["plt_yRange"]     = [0.0,40.0]
    config["plt_GaussFilt"]  = 3.0
    config["xMajor_Nticks"]  = 5
    config["yMajor_Nticks"]  = 5
    config["MinimalOut"]     = True
    config["plt_LegFontSize"]= 13
    config["plt_LegLocation"]= "lower=right"

    # ---------------------------------------- #
    # --- [4] プロット                     --- #
    # ---------------------------------------- #
    fig     = pl1.plot1D( FigName=pngFile, config=config )
    for hData in Data:
        fig.addPlot( xAxis=hData["ptime"], yAxis=hData["x2Xpt"], label=Legs[ hData["job"] ] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    
# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    pngDir  = "./png/rXptPlot/"
    jobBase = None
    ret     = rXptPlotter( jobBase=jobBase, pngDir=pngDir )
