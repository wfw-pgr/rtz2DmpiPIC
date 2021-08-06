import numpy                  as np
import myStyle.LoadConfig     as lcf
import myStyle.StackBarPlot   as sbp
import myStyle.configSettings as cfs
import myConvert.arr2dct      as a2d

# ======================================== #
# ===  棒グラフメーカー                === #
# ======================================== #
def singleBarMaker( config=None, Flag__Display=True ):
    # ---------------------------------------- #
    # --- [1]    引数チェック              --- #
    # ---------------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    # ---------------------------------------- #
    # --- [2] ユーザ指定領域               --- #
    # ---------------------------------------- #
    job     = "Init21_"
    datFile = "{0}/EnergyBar.dat".format( job )
    pngFile = "{0}/EnergyBar.png".format( job )
    pkeys   = ["UMFx","UMFy","UMFz","UEFx","UEFy","UEFz","Ukex","Ukey","Ukez","Ukix","Ukiy","Ukiz",]
    colors  = ["green","paleGreen","darkCyan","Gold","Orange","yellow",\
               "RoyalBlue","mediumBlue","Cyan","magenta","Crimson","purple" ]
    # ---------------------------------------- #
    # --- [3] データ読み取り               --- #
    # ---------------------------------------- #
    Data  = ( np.loadtxt( datFile, comments="#" ) ).reshape( (-1,1))
    with open( datFile, "r" ) as f:
        items = ( ( f.readline() ).replace("#","") ).split()
    Data  = a2d.arr2dct( Data=Data, keys=items )
    pData = ( np.array( [ Data[key] for key in pkeys ] ) ).reshape( (1,-1))
    coef  = 100.0 / np.sum( pData )
    if ( Flag__Display ):
        for key in pkeys: print( "[singleBarMaker] {0} :: {1}  ( {2} % )".format( key, Data[key], Data[key]*coef ) )
    # ---------------------------------------- #
    # --- [2] コンフィグ & プロット        --- #
    # ---------------------------------------- #
    cfs.configSettings( config=config, configType="cMap2D_thesis" )
    config["MinimalOut"] = True
    config["FigSize"]    = (3,10)
    sbp.StackBarPlot( Data   =pData  , BarOccupationRatio=1.0, colors=colors, \
                      FigName=pngFile, config=config )

# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    singleBarMaker()
