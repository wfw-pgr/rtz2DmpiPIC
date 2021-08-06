import numpy                   as np
import myStyle.LoadConfig      as lcf
import mPICsee.TimeUnitConvert as tuc
import myStyle.plot1D          as pl1
import myStyle.configSettings  as cfs
import myConvert.arr2dct       as a2d

# ======================================== #
# ===  合体率 ( mRate )  を 図示       === #
# ======================================== #
def mRatePlotter( jobBase=None, pngDir=None, pngFile=None, config=None ):
    # ------------------------------- #
    # --- [1]   引数チェック      --- #
    # ------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( jobBase is None ): jobBase = config["pic_jobDir"]
    if ( pngDir  is None ): pngDir  = "./"
    if ( pngFile is None ): pngFile = "{0}mRatePlot.png".format( pngDir )
    # ------------------ #
    # -- 対象ジョブ名 -- #
    # ------------------ #
    jobs  = ["CtrO21_01", "CtrI21_01", "Co-21_01"]
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
    config["plt_xAutoRange"] = False
    config["plt_xRange"]     = [0.0,120.0]
    config["plt_yAutoRange"] = False
    config["plt_yRange"]     = [0.0,1.2]
    config["plt_GaussFilt"]  = 2.0
    config["xTitle"]         = ""
    config["yTitle"]         = ""
    config["xMajor_Nticks"]  = 5
    config["yMajor_Nticks"]  = 7
    config["MinimalOut"]     = True
    config["cursor_y"]       = 1.0
    config["cursor_style"]   = "-"
    config["plt_LegFontSize"]= 13
    config["plt_LegLocation"]= "lower=right"

    # ---------------------------------------- #
    # --- [4] プロット                     --- #
    # ---------------------------------------- #
    fig     = pl1.plot1D( FigName=pngFile, config=config )
    for hData in Data:
        fig.addPlot( xAxis=hData["ptime"], yAxis=hData["mergR"], label=Legs[ hData["job"] ] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()

    
# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    pngDir  = "./png/mRatePlot/"
    jobBase = None
    ret     = mRatePlotter( jobBase=jobBase, pngDir=pngDir )
