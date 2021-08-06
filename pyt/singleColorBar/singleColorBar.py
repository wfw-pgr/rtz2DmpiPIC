import sys
import numpy              as np
import myStyle.LoadConfig as lcf
import myStyle.cMap2D     as clm

def singleColorBar( MaxMin    =None , Nlevels =None, FigSize=None, OutFile     =None , \
                    horizontal=None , vertical=None, Aspect =None, AdjustAspect=1.5  , \
                    MinimalOut=True , config  =None  ):
    # ------------------------ #
    # --- [1] 引数チェック --- #
    # ------------------------ #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( MaxMin  is None ): MaxMin  = config["cmp_MaxMin"]
    if ( Nlevels is None ): Nlevels = config["cmp_Nlevels"]
    if ( OutFile is None ): OutFile = "ColorBar.png"
    if ( FigSize is None ):
        if   ( horizontal  ):
            FigSize  = (5,1)
        elif ( vertical    ):
            FigSize  = (1,5)
        else:
            vertical, FigSize = True, (1,5)
    print( FigSize )
    # ------------------------ #
    # --- [2] バー作成     --- #
    # ------------------------ #
    clb  = np.zeros( (2,Nlevels) )
    clbs = np.linspace( MaxMin[0], MaxMin[1], Nlevels )
    clb[0,:], clb[1,:] = clbs, clbs
    
    # --------------------------- #
    # --- [3] 設定 & プロット --- #
    # --------------------------- #
    if ( Aspect is None ): Aspect = float( FigSize[1] ) / float( FigSize[0] )
    config["xTitle"]     = ""
    config["yTitle"]     = ""
    config["clb_sw"]     = False
    config["FigSize"]    = FigSize
    config["FigName"]    = OutFile
    config["MinimalOut"] = MinimalOut
    config["xMajor_NoLabel"] = True
    config["yMajor_NoLabel"] = True
    config["cmp_ColorTable"] = "hot"
    if   ( horizontal ):
        config["axes_y_Nticks"] = 0
        config["cmp_Position" ] = [0.15,0.45,0.90,0.90]
        widt = np.linspace( 0.0, ( MaxMin[1]-MaxMin[0] )*(Aspect/AdjustAspect), 2 )
        clm.cMap2D( xAxis =clbs   , yAxis =widt, cMap=np.transpose( clb ), config=config )
    elif ( vertical   ):
        config["axes_x_Nticks"] = 0
        config["cmp_Position" ] = [0.45,0.15,0.90,0.90]
        widt = np.linspace( 0.0, ( MaxMin[1]-MaxMin[0] )/(Aspect*AdjustAspect), 2 )
        clm.cMap2D( xAxis =widt   , yAxis =clbs, cMap=clb, config=config )
    else:
        sys.exit( "[singleColorBar] vertical=False & horizontal=False" )


# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    MaxMin  = args["MaxMin"]
    horizontal, vertical = None, None
    FigSize = (8,1)
    if ( args["key"][0] == "horizontal" ): horizontal, vertical=True , False
    if ( args["key"][0] == "vertical"   ): horizontal, vertical=False, True
    print( "horizontal :: {0}".format( horizontal ) )
    print( "vertical   :: {0}".format( vertical   ) )
    print( "FigSize    :: {0}".format( FigSize    ) )
    singleColorBar( MaxMin=MaxMin, Nlevels=10, FigSize=FigSize, \
                    horizontal=horizontal, vertical=vertical    )
