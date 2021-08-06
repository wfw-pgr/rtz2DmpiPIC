import numpy                  as np
import myStyle.plot1D         as pl1
import mPICsee.FetchPIC       as fpc
import myStyle.LoadConfig     as lcf
import myStyle.configSettings as cfs

# --- pic --- #
def pic2slc( job    =None, kstep  =None, jobDir   =None , OutFile=None, OutDir =None, \
             keys   =None, labels =None, pltAxis  =None , x1Range=None, x2Range=None, writeData=False, \
             config =None  ):
    # -------------------------- #
    # --- [1] 引数チェック   --- #
    # -------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}png/".format( jobDir )
    if ( OutFile is None ): OutFile = "Field1D{0:08}".format( kstep )
    if ( pltAxis is None ): pltAxis = "x2"
    if ( keys    is None ): keys    = [config["pic_key"]]
    if ( labels  is None ): labels  = [None]*len( keys )

    # -------------------------- #
    # --- [2] データ呼び出し --- #
    # -------------------------- #
    #  -- 軸の選択 / 平均化  --  #
    if   ( pltAxis == "x1" ):
        meanAxis = 0
        if ( x1Range is None ): x1Range = config["Data_x2Range" ]
        if ( x2Range is None ): x2Range = config["Data_avgRange"]
    elif ( pltAxis == "x2" ):
        meanAxis = 1
        if ( x1Range is None ): x1Range = config["Data_avgRange"]
        if ( x2Range is None ): x2Range = config["Data_x1Range" ]
    else:
        import sys
        sys.exit( "[pic2slc]  ERROR!! :: plt_Axis??" )
    #  -- データ読み出し -- #
    config["pic_Temperature"] = True
    keys = ["pexx", "peyy", "pezz"]
    Data = fpc.FetchPIC( job    =job    , jobDir =jobDir , kstep =kstep, \
                         keys   =keys   , config =config  )
    #, \x1Range=x2Range, x2Range=x1Range )
    Tprof = ( Data["pexx"] + Data["peyy"] + Data["pezz"] )/ 3.0
    print( Tprof.shape )

    # ----------------------- #
    # --- [3] プロット    --- #
    # ----------------------- #
    #  -- [3-1] プロット  --  #
    iXpt = 512
    jXpt = 330
    print( Data["xAxis"][iXpt] )
    print( Data["yAxis"][jXpt] )
    cfs.configSettings( configType="plot1D_def", config=config )
    xAxis1   = Data["xAxis"]
    xAxis2   = Data["yAxis"]
    yAxis1   = np.ravel( np.mean( Tprof[(jXpt-7):(jXpt+7),:], axis=0 ) )
    yAxis2   = np.ravel( Tprof[:,iXpt])
    pl1.plot1D( xAxis=xAxis1, yAxis=yAxis1, FigName="Te1.png", config=config )
    pl1.plot1D( xAxis=xAxis2, yAxis=yAxis2, FigName="Te2.png", config=config )


# --------------------------------------- #
if ( __name__=="__main__" ):
    # -- 引数 -- #
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    OutDir  = "./"
    config  = lcf.LoadConfig()
    job     = args["job"]
    kstep   = args["kstep"]
    jobDir  = args["jobDir"]
    keys    = args["key"]

    # --- コンフィグの設定 --- #
    config["FigSize"]         = (6,2)
    config["plt_yAutoRange"]  = False
    config["plt_yRange"]      = [0.0,0.15]
    config["plt_color"]       = "indigo"
    config["AutoRange"]       = False
    config["AutoTicks"]       = True
    config["axes_x_Nticks"]   = 5
    config["axes_y_Nticks"]   = 5
    config["MinimalOut"]      = True
    config["axes_x_off"]      = True
    config["axes_y_off"]      = True
    # -- 描く -- #
    pltAxis = "x1"
    x1Range = [-8.6,+8.6]
    x2Range = [ 5.0,+5.5]
    print( "[pic2slc] ( job, kstep ) = ( {0}, {1} ) is under processing...".format( job, kstep ) )
    vtk     = pic2slc( job    =job    , jobDir=jobDir, kstep =kstep , \
                       OutDir =OutDir , keys  =keys  , config=config, \
                       pltAxis=pltAxis, x1Range=x1Range, x2Range=x2Range )
    print( "[pic2slc]   Done" )
