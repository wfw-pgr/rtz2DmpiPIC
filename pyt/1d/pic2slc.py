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
    if ( pltAxis is None ): pltAxis = config["plt_slcAxis"]
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
    Data         = fpc.FetchPIC( job    =job    , jobDir =jobDir , kstep =kstep, \
                                 keys   =keys   , config =config , \
                                 x1Range=x1Range, x2Range=x2Range )

    # ----------------------- #
    # --- [3] プロット    --- #
    # ----------------------- #
    #  -- [3-1] プロット  --  #
    cfs.configSettings( configType="plot1D_def", config=config )
    fig          = pl1.plot1D( FigName=OutFile, config=config )
    for i, key in enumerate(keys):
        xAxis    = Data["xAxis"] if ( pltAxis == "x1" ) else Data["yAxis"]
        yAxis    = np.ravel( np.mean( Data[key], axis=meanAxis ) )
        fig.addPlot( xAxis=xAxis, yAxis=yAxis, label=labels[i] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    #  -- [3-2] .dat 書き出し --  #
    if ( writeData ):
        datFile = OutFile.replace( ".png", ".dat" )
        outData = np.concatenate( [xAxis[np.newaxis,:], yAxis[np.newaxis,:]], axis=0 )
        np.savetxt( datFile, outData )
        print( "[pic2slc] {0} is outputed...".format( datFile ) )


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
    config["FigSize"]         = (6,3)
    config["AutoRange"]       = True
    config["AutoTicks"]       = True
    config["cmp_xRange"]      = [-0.5,+0.5]
    config["cmp_yRange"]      = [-0.5,+0.5]
    config["axes_x_Nticks"]   = 5
    config["axes_y_Nticks"]   = 5
    config["cnt_Nlevels"]     = 50
    config["cnt_thick"]       = 0.10
    config["cmp_AutoLevel"]   = True
    config["cmp_MaxMin"]      = [0.0,0.0020]
    config["clb_sw"]          = False
    config["MinimalOut"]      = False
    # -- 描く -- #
    print( "[pic2slc] ( job, kstep ) = ( {0}, {1} ) is under processing...".format( job, kstep ) )
    vtk     = pic2slc( job    =job    , jobDir=jobDir, kstep =kstep , \
                       OutDir =OutDir , keys  =keys  , config=config )
    print( "[pic2slc]   Done" )
