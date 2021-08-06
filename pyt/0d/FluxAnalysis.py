import numpy                         as np
import myStyle.LoadConfig            as lcf
import myUtils.progressBar           as pgb
import mPICsee.FetchPIC              as fpc
import myConvert.arr2dct             as a2d
import myStyle.plot1D                as pl1
import myStyle.configSettings        as cfs
import mPICsee.TimeUnitConvert       as tuc
import myAnalysis.createMask         as msk
import myBasicAlgs.robustInv         as inv
import myBasicAlgs.cylindricalVolume as cvl
import myStyle.pickColor             as pcl

# =========================================================== #
# ===  フラックス / Ip 解析用のルーチン ( ハンドラ )      === #
# =========================================================== #
def FluxAnalysis( job    =None, jobDir=None,  ksteps=None , retSize=21  , \
                  OutFile=None, OutDir=None,  config=None , pngDir =None, CnsFile=None, \
                  NewFile=True, Flag__CallCalculater=False, Flag__CallPlotter=True ):
    # ------------------------------- #
    # --- [1]   引数チェック      --- #
    # ------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( ksteps  is None ): ksteps  = np.arange( config["arg_Iter"] )*config["arg_Step"]+config["arg_Init"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}dat/".format( jobDir )
    if ( pngDir  is None ): pngDir  = "{0}png/".format( jobDir )
    if ( OutFile is None ): OutFile = OutDir + "FluxAnalysis.dat"
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    # ------------------------------- #
    # --- [2]  全磁束 / 全電流    --- #
    # ------------------------------- #
    if ( Flag__CallCalculater ):
        Nsteps = len( ksteps )
        if ( NewFile ): ( open( OutFile, "w" ) ).close()    # -- 新規ファイル -- #
        pgb.progressBar( iteration=0, total=Nsteps )        # -- 進捗表示     -- #
        for i,kstep in enumerate( ksteps ):
            ret = FluxAnalysis_calculater( job    =job    , jobDir=jobDir, kstep=kstep, \
                                           config=config )
            pgb.progressBar( iteration=i+1, total=Nsteps )
            with open( OutFile, "ab" ) as f:
                np.savetxt( f, np.array( ret ).reshape( (1,retSize) ) )
    # ------------------------------- #
    # --- [3] 返却 /  描画        --- #
    # ------------------------------- #
    #  -- [3-1] .png 書き出し     --  #
    if ( Flag__CallPlotter ):
        FluxAnalysis_plotter( datFile=OutFile, config=config, pngDir=pngDir, CnsFile=CnsFile )


def FluxAnalysis_calculater( job=None, jobDir=None, kstep=None, config=None ):
    # -------------------------- #
    # --- [1] データ呼び出し --- #
    # -------------------------- #
    #  -- [1-1] time / rAxis --  #
    fstep   = float( kstep )
    keys    = ["psi","By","Jy","ne","ni"]
    Data    = fpc.FetchPIC( job=job, jobDir=jobDir, kstep=kstep, keys=keys, config=config )
    mask    = msk.createMask( Flux=Data["psi"], Flag__privateFluxMode=True )
    area    = ( Data["xAxis"][1]-Data["xAxis"][0] )*( Data["yAxis"][1]-Data["yAxis"][0] )
    volm    = cvl.cylindricalVolume( zAxis=Data["xAxis"], rAxis=Data["yAxis"] )
    #  -- [1-3] マスク面積分 --  #
    ByPrv0  = np.sum( Data["By"] * mask["fMask0"] * area )
    ByPrv1  = np.sum( Data["By"] * mask["fMask1"] * area )
    ByPrv2  = np.sum( Data["By"] * mask["fMask2"] * area )
    ByComn  = np.sum( Data["By"] * mask["fMask3"] * area )
    JyPrv0  = np.sum( Data["Jy"] * mask["fMask0"] * area )
    JyPrv1  = np.sum( Data["Jy"] * mask["fMask1"] * area )
    JyPrv2  = np.sum( Data["Jy"] * mask["fMask2"] * area )
    JyComn  = np.sum( Data["Jy"] * mask["fMask3"] * area )
    nePrv0  = np.sum( Data["ne"] * mask["fMask0"] * volm )
    niPrv0  = np.sum( Data["ni"] * mask["fMask0"] * volm )
    neComn  = np.sum( Data["ne"] * mask["fMask3"] * volm )
    niComn  = np.sum( Data["ni"] * mask["fMask3"] * volm )
    ByTot   = np.sum( Data["By"] * area )
    JyTot   = np.sum( Data["Jy"] * area )
    neTot   = np.sum( Data["ne"] * volm )
    niTot   = np.sum( Data["ni"] * volm )
    csArea  = np.sum( area * mask["fMask3"] )
    volume  = np.sum( volm * mask["fMask3"] )
    psiMax  = Data["psiMax"]
    time    = Data["time"]
    return( [ fstep , time  , ByPrv0, ByPrv1, ByPrv2, ByComn, JyPrv0, JyPrv1, JyPrv2, JyComn, \
              nePrv0, niPrv0, neComn, niComn, ByTot , JyTot , neTot , niTot , csArea, volume, psiMax ] )

# --- .png 書き出し --- #
def FluxAnalysis_plotter( datFile=None, config=None, pngDir=None, CnsFile=None ):
    # ---------------------------------------- #
    # --- [1] プロット準備                  -- #
    # ---------------------------------------- #
    #  -- [1-1] データ読込                 --  #
    Data  = np.transpose( np.loadtxt( datFile ) )
    keys  = ["fstep" , "time"  , "ByPrv0", "ByPrv1", "ByPrv2", "ByComn", \
             "JyPrv0", "JyPrv1", "JyPrv2", "JyComn", "nePrv0", "niPrv0", "neComn", "niComn", \
             "ByTot" , "JyTot" , "neTot" , "niTot" , "csArea", "volume", "psiMax" ]
    Data  = a2d.arr2dct( Data=Data, keys=keys )
    ptime = tuc.TimeUnitConvert( Data=Data["time"], unit="wpi", CnsFile=CnsFile )
    #  -- [1-1] コンフィグ                 --  #
    cfs.configSettings( config=config, configType="plot1D_lateral" )
    cfs.configSettings( config=config, configType="plot1D_stack"   )
    cfs.configSettings( config=config, configType="FilterClear"    )
    config["plt_MedianFilt"]  = 3
    config["plt_GaussFilt" ]  = 1.2
    config["plt_xAutoRange"]  = False
    config["plt_xRange"]      = [0.,120.]
    config["plt_yAutoRange"]  = False
    config["plt_LegFontSize"] = 14
    colors  = pcl.pickColor( nColor=20, pallete="mild", shuffle=False )
    labels  = {"ByPrv1":"$\Phi_{sph1}$", "ByPrv2":"$\Phi_{sph2}$","ByComn":"$\Phi_{com}$","psiMax":"$\Psi_{Max}$", \
               "JyPrv1":"$I_{sph1}$", "JyPrv2":"$I_{sph2}$", "JyComn":"$I_{p}$", \
               "neDiff":"$n_{e}^{Down}$", "nePrv0":"$n_{e}^{Up}$", "neComn":"$n_{e}^{com}$" }
    
    # ---------------------------------------- #
    # --- [2] プロット                     --- #
    # ---------------------------------------- #
    #  --------------------  #
    #  --  Toroidal Flux --  #
    #  --------------------  #
    OutFile = pngDir + "t-ToroidalFlux.png"
    psiMaxs = np.abs( Data["psiMax"] )
    config["plt_yRange"] = [-1000.0,1000.0]
    config["plt_LegLocation"] = "lower right"
    config["plt_LegNColumn"]  = 2
    fig = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=ptime, yAxis=Data["ByPrv1"], color=colors[1], label=labels["ByPrv1"] )
    fig.addPlot( xAxis=ptime, yAxis=Data["ByPrv2"], color=colors[0], label=labels["ByPrv2"] )
    fig.addPlot( xAxis=ptime, yAxis=Data["ByComn"], color=colors[3], label=labels["ByComn"] )
    fig.addPlot( xAxis=ptime, yAxis=psiMaxs       , color=colors[4], label=labels["psiMax"] )
    fig.setAxis()
    fig.addLegend()
    fig.writeFigure()
    #  --------------------  #
    #  -- plasma Current --  #
    #  --------------------  #
    OutFile = pngDir + "t-PlasmaCurrent.png"
    config["plt_yRange"] = [0.0,100.0]
    config["plt_LegLocation"] = "upper right"
    config["plt_LegNColumn"]  = 3
    fig = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=ptime, yAxis=Data["JyPrv1"], color=colors[2] , label=labels["JyPrv1"] )
    fig.addPlot( xAxis=ptime, yAxis=Data["JyPrv2"], color=colors[15], label=labels["JyPrv2"] )
    fig.addPlot( xAxis=ptime, yAxis=Data["JyComn"], color=colors[4] , label=labels["JyComn"] )
    fig.setAxis()
    fig.addLegend()
    fig.writeFigure()
    #  --------------------  #
    #  --    ne Total    --  #
    #  --------------------  #
    config["plt_yRange"]      = [0.0,0.6]
    config["plt_LegLocation"] = "upper right"
    config["plt_LegNColumn"]  = 3
    nePrv0  = Data["nePrv0"] * inv.robustInv( Data["neTot"] )
    neComn  = Data["neComn"] * inv.robustInv( Data["neTot"] )
    neDiff  = neComn - nePrv0
    OutFile = pngDir + "t-neTotal.png"
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=ptime, yAxis=neDiff, color=colors[7] , label=labels["neDiff"]  )
    fig.addPlot( xAxis=ptime, yAxis=nePrv0, color=colors[8] , label=labels["nePrv0"]  )
    fig.addPlot( xAxis=ptime, yAxis=neComn, color=colors[10], label=labels["neComn"]  )
    fig.setAxis()
    fig.addLegend()
    fig.writeFigure()
    #  --------------------  #
    #  --  Volume / Area --  #
    #  --------------------  #
    config["plt_yAutoRange"] = True
    OutFile = pngDir + "t-VolumeArea.png"
    pl1.plot1D( xAxis=ptime, yAxis=Data["volume"], config=config, FigName=OutFile )
    OutFile = pngDir + "t-csArea.png"
    pl1.plot1D( xAxis=ptime, yAxis=Data["csArea"], config=config, FigName=OutFile )
    #  --------------------  #
    #  --  Volume / Area --  #
    #  --------------------  #
    config["plt_yAutoRange"] = False
    config["plt_yRange"]     = [500,700]
    OutFile = pngDir + "t-psiMax.png"
    pl1.plot1D( xAxis=ptime, yAxis=psiMaxs, config=config, FigName=OutFile )


# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    pngDir  = "./png/FluxAnalysis/"
    print( "[FluxAnalysis] {0} is under processed...".format( args["job"] ), end="\n" )
    Output  = FluxAnalysis( job=args["job"], pngDir=pngDir, ksteps=args["ksteps"] )
    print( "\t [Done]" )






# rAxis   = np.linspace(   const["x2Min"]     , const["x2Min"]+const["x2Leng"], const["LJ"] )
# zAxis   = np.linspace( - const["x1Leng"]*0.5, const["x1Leng"]*0.5           , const["LI"] )
#  -- [1-2] データ取得   --  #
# Data    = m2o.mpi2one( job=job, jobDir=jobDir, kstep=kstep, config=config )
# keyDict = {"By":1,"Bz":2,"Jy":3,"ne":18,"ni":19}
# Bz      = np.ascontiguousarray( Data[ keyDict["Bz"],:,:] )
# psiMax  = [1.0]
# Flux    = cfl.calcFlux ( rg=rAxis, Bz=Bz , indexing='ji', Flag__normalize=True, fluxMax=psiMax )
# import fLIB.fLIB__mpi2one      as m2o
# import mPICsee.calcFlux        as cfl
