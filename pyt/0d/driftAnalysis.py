import numpy                         as np
import mPICsee.FetchPIC              as fpc
import mPICsee.driftField            as dfd
import mPICsee.TimeUnitConvert       as tuc
import myStyle.plot1D                as pl1
import myStyle.LoadConfig            as lcf
import myStyle.configSettings        as cfs
import myUtils.progressBar           as pgb
import myBasicAlgs.cylindricalVolume as cvl
import myConvert.arr2dct             as a2d
import myAnalysis.outletMask         as olm

# ======================================== #
# ===    Jd.Eの時間領域 解析           === #
# ======================================== #
def driftAnalysis( job    =None, jobDir=None, ksteps =None , \
                   datFile=None, datDir=None, pngFile=None , pngDir =None, config=None, \
                   NewFile=True, Flag__CallCalculater=True , Flag__CallPlotter=True, ptype="e", RawData=True ):
    # ------------------------------- #
    # --- [1]   引数チェック      --- #
    # ------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( ksteps  is None ): ksteps  = np.arange( config["arg_Iter"] )*config["arg_Step"]+config["arg_Init"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( datDir  is None ): datDir  = "{0}dat/".format( jobDir )
    if ( pngDir  is None ): pngDir  = "{0}png/".format( jobDir )
    if ( datFile is None ): datFile = "{0}EJ{1}_driftAnalysis.dat".format( datDir, ptype )
    if ( pngFile is None ): pngFile = "{0}EJ{1}_driftAnalysis.png".format( pngDir, ptype )
    # ------------------------------- #
    # --- [2]  J.E を解析         --- #
    # ------------------------------- #
    if ( Flag__CallCalculater ):
        Nsteps = len( ksteps )
        if ( NewFile ):                                # -- 新規ファイル -- #
            items = "# fstep ptime sEJ_ExB1 sEJ_grB1 sEJ_cvt1 sEJ_plr1 sEJ_ExB2 sEJ_grB2 sEJ_cvt2 sEJ_plr2" + "\n"
            with open( datFile, "w" ) as f: f.write( items )
        pgb.progressBar( iteration=0, total=Nsteps )   # -- 進捗表示     -- #
        for i,kstep in enumerate( ksteps ):
            ret = EdJ_calculator( job=job, jobDir=jobDir, kstep=kstep, config=config, ptype=ptype, RawData=RawData )
            pgb.progressBar( iteration=i+1, total=Nsteps )
            with open( datFile, "ab" ) as f:
                np.savetxt( f, np.array( ret ).reshape( (1,-1) ) )
    # ------------------------------- #
    # --- [3] 返却 /  描画        --- #
    # ------------------------------- #
    #  -- [3-1] .png 書き出し     --  #
    if ( Flag__CallPlotter ):
        driftAnalysis_plotter( datFile=datFile, config=config, pngFile=pngFile )


# ======================================== #
# ===       Jd.Eを計算して表示         === #
# ======================================== #
def EdJ_calculator( job  =None, jobDir=None, kstep =None, config=None, ptype="e", RawData=True ):
    # -------------------------- #
    # --- [1]  引数チェック  --- #
    # -------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"    .format( config["pic_jobDir"], job )
    # -------------------------- #
    # --- [2] データ呼び出し --- #
    # -------------------------- #
    #  -- [2-1] 読込         --  #
    config["cmp_nFilter"]    = 10
    config["cmp_LinearFilt"] = 0.2
    fstep = float( kstep )
    keys  = [ "Ex", "Ey", "Ez", "ne", "ni", "psi" ]
    Data  = fpc.FetchPIC( job =job , jobDir =jobDir , kstep =kstep , FilterKeys=keys, \
                          keys=keys, RawData=RawData, config=config )
    xAxis = Data["xAxis"]
    yAxis = Data["yAxis"]
    flux  = Data["psi"  ]
    drift = dfd.driftField     ( job=job, jobDir=jobDir, kstep=kstep, RawData=RawData, ptype=ptype, config=config )
    ptime = tuc.TimeUnitConvert( job=job, jobDir=jobDir, kstep=kstep, unit="wpi" )
    #  -- [2-2] 準備           --  #
    maskDict = olm.outletMask( Flux=flux )
    mask1    = maskDict["outlet1"]
    mask2    = maskDict["outlet2"]
    volm     = cvl.cylindricalVolume( zAxis=xAxis, rAxis=yAxis )
    # ---------------------------- #
    # --- [3]  E.Jd            --- #
    # ---------------------------- #
    if   ( ptype=="e" ):
        qns  = - 1.0 * Data["ne"]
    elif ( ptype=="i" ):
        qns  = + 1.0 * Data["ni"]
    qnvolm1  = qns * mask1 * volm
    qnvolm2  = qns * mask2 * volm
    sEJ_ExB1 = np.sum( ( drift["vExBx"]*Data["Ex"] + drift["vExBy"]*Data["Ey"] + drift["vExBz"]*Data["Ez"] ) * qnvolm1 )
    sEJ_grB1 = np.sum( ( drift["vgrBx"]*Data["Ex"] + drift["vgrBy"]*Data["Ey"] + drift["vgrBz"]*Data["Ez"] ) * qnvolm1 )
    sEJ_cvt1 = np.sum( ( drift["vcvtx"]*Data["Ex"] + drift["vcvty"]*Data["Ey"] + drift["vcvtz"]*Data["Ez"] ) * qnvolm1 )
    sEJ_plr1 = np.sum( ( drift["vplrx"]*Data["Ex"] + drift["vplry"]*Data["Ey"] + drift["vplrz"]*Data["Ez"] ) * qnvolm1 )
    sEJ_ExB2 = np.sum( ( drift["vExBx"]*Data["Ex"] + drift["vExBy"]*Data["Ey"] + drift["vExBz"]*Data["Ez"] ) * qnvolm2 )
    sEJ_grB2 = np.sum( ( drift["vgrBx"]*Data["Ex"] + drift["vgrBy"]*Data["Ey"] + drift["vgrBz"]*Data["Ez"] ) * qnvolm2 )
    sEJ_cvt2 = np.sum( ( drift["vcvtx"]*Data["Ex"] + drift["vcvty"]*Data["Ey"] + drift["vcvtz"]*Data["Ez"] ) * qnvolm2 )
    sEJ_plr2 = np.sum( ( drift["vplrx"]*Data["Ex"] + drift["vplry"]*Data["Ey"] + drift["vplrz"]*Data["Ez"] ) * qnvolm2 )
    # ---------------------------- #
    # --- [5]     返却         --- #
    # ---------------------------- #
    return( [ fstep, ptime, sEJ_ExB1, sEJ_grB1, sEJ_cvt1, sEJ_plr1, sEJ_ExB2, sEJ_grB2, sEJ_cvt2, sEJ_plr2 ] )


def driftAnalysis_plotter( datFile=None, pngFile=None, config=None ):
    # ------------------------- #
    # --- [1] 準備          --- #
    # ------------------------- #
    Data = np.transpose( np.loadtxt( datFile, comments="#" ) )
    with open( datFile, "r" ) as f:
        items = ( ( f.readline() ).replace("#","") ).split()
    Legp = [ r"$E \times B$", r"$\nabla B$", "curvature", "polarization" ]
    Legs = [ "fstep","Time"] + Legp + Legp
    Data = a2d.arr2dct( Data=Data, keys=items )
    Legs = a2d.arr2dct( Data=Legs, keys=items )
    cfs.configSettings( config=config, configType="plot1D_lateral" )
    cfs.configSettings( config=config, configType="FilterClear"    )
    # ------------------------- #
    # --- [2] プロット      --- #
    # ------------------------- #
    #  --   E.Jd   --  #
    fig     = pl1.plot1D( FigName=pngFile, config=config )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sEJ_ExB1"], label=Legs["sEJ_ExB1"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sEJ_grB1"], label=Legs["sEJ_grB1"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sEJ_cvt1"], label=Legs["sEJ_cvt1"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sEJ_plr1"], label=Legs["sEJ_plr1"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sEJ_ExB2"], label=Legs["sEJ_ExB2"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sEJ_grB2"], label=Legs["sEJ_grB2"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sEJ_cvt2"], label=Legs["sEJ_cvt2"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sEJ_plr2"], label=Legs["sEJ_plr2"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()

# ======================================== #
# ===           テスト実行             === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    OutDir  = "./png/driftAnalysis/"
    driftAnalysis( job   =args["job"]   , jobDir=args["jobDir"], \
                   ksteps=args["ksteps"] )
