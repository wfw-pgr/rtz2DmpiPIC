import numpy                  as np
import myStyle.LoadConfig     as lcf
import myStyle.plot1D         as pl1
import myStyle.configSettings as cfs
import myConvert.arr2dct      as a2d
import myConvert.pilearr      as pil

# ======================================== #
# ===  プロッター                      === #
# ======================================== #
def mrUpt_plotter( config=None, pngDir=None, pngFile=None, Flag__process="delta" ):
    # ---------------------------------------- #
    # --- [1] 引数チェック                 --- #
    # ---------------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( pngDir  is None ): pngDir  = "./png/"
    if ( pngFile is None ): pngFile = "mrUpt_{0}.png".format( Flag__process )
    # ---------------------------------------- #
    # --- [2] データ読込                   --- #
    # ---------------------------------------- #
    #  --  読込  -- #
    InpFile1 = "mr-Upt_CtrO.dat"
    InpFile2 = "mr-Upt_Co.dat"
    Data1    = np.transpose( np.loadtxt( InpFile1, comments="#" ) )
    Data2    = np.transpose( np.loadtxt( InpFile2, comments="#" ) )
    #  --  CtrO  -- #
    idx1     = [[0,2,4], [1,3,5]]
    Data1b   = Data1[ 3:, idx1[0] ]
    Data1a   = Data1[ 3:, idx1[1] ]
    #  --  CoMg  -- #
    idx2     = [[0,2],[1,3]]
    Data2b   = Data2[ 3:, idx2[0] ]
    Data2a   = Data2[ 3:, idx2[1] ]
    if ( Flag__process == "ratio" ):
        Data1d  = Data1a / Data1b
        Data2d  = Data2a / Data2b
    if ( Flag__process == "delta" ):
        Data1d  = ( Data1a - Data1b ) / Data1b
        Data2d  = ( Data2a - Data2b ) / Data2b
    Data1r   = pil.pilearr( [Data1[0:3,idx1[1]], Data1d], NoNewAxis=True, axis=0 )
    Data2r   = pil.pilearr( [Data2[0:3,idx2[1]], Data2d], NoNewAxis=True, axis=0 )
    # --- プロット用データまとめ --- #
    keys     = ["kstep", "Mr", "volm", "ethe", "ethi", "eFle", "eFli", "epte", "epti"]
    CtrO     = a2d.arr2dct( Data=Data1r, keys=keys )
    CoMg     = a2d.arr2dct( Data=Data2r, keys=keys )
    if   ( Flag__process=="ratio" ):
        labels   = {"CtrO_epte":"$U_{e}^{Ctr}$", "CtrO_epti":"$U_{i}^{Ctr}$", \
                    "CoMg_epte":"$U_{e}^{Co}$" , "CoMg_epti":"$U_{i}^{Co}$"   }
    elif ( Flag__process=="delta" ):
        labels   = {"CtrO_epte":"$\Delta U_{e}^{Ctr}$", "CtrO_epti":"$\Delta U_{i}^{Ctr}$", \
                    "CoMg_epte":"$\Delta U_{e}^{Co}$" , "CoMg_epti":"$\Delta U_{i}^{Co}$"   }
    # ---------------------------------------- #
    # --- [3] プロット                     --- #
    # ---------------------------------------- #
    # -- コンフィグ -- #
    cfs.configSettings( config=config, configType="plot1D_def" )
    config["plt_yAutoRange"] = False
    config["plt_xAutoRange"] = False
    config["plt_linewidth"]  = 1.0
    config["plt_marker"]     = "D"
    config["plt_xRange"]     = [0.0,500.0]
    config["plt_yRange"]     = [0.0,12.0 ]
    config["xMajor_Nticks"]  = 6
    config["xTitle"]         = ""
    config["yTitle"]         = ""
    config["plt_LegNColumn"] = 2
    # --  プロット  -- #
    fig      = pl1.plot1D( config=config, FigName=pngFile )
    fig.addPlot( xAxis=CtrO["Mr"], yAxis=CtrO["epte"], label=labels["CtrO_epte"] )
    fig.addPlot( xAxis=CtrO["Mr"], yAxis=CtrO["epti"], label=labels["CtrO_epti"] )
    fig.addPlot( xAxis=CoMg["Mr"], yAxis=CoMg["epte"], label=labels["CoMg_epte"] )
    fig.addPlot( xAxis=CoMg["Mr"], yAxis=CoMg["epti"], label=labels["CoMg_epti"] )
    fig.setAxis()
    fig.addLegend()
    fig.writeFigure()
    
# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    pngDir = "./png/"
    ret    = mrUpt_plotter( pngDir=pngDir )
