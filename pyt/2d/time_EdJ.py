import numpy                         as np
import myStyle.plot1D                as pl1
import mPICsee.FetchPIC              as fpc
import myStyle.LoadConfig            as lcf
import myBasicAlgs.robustInv         as riv
import myStyle.configSettings        as cfs
import myBasicAlgs.cylindricalVolume as cvl

# ======================================== #
# ===    J.Eの時間領域 解析            === #
# ======================================== #
def EdJAnalysis( job    =None, jobDir=None, kstep  =None , CnsFile=None, uorJ  =None, \
                 datFile=None, datDir=None, pngFile=None , pngDir =None, config=None, \
                 NewFile=True, Flag__CallCalculater=False, Flag__CallPlotter=True, ):
    # ------------------------------- #
    # --- [1]   引数チェック      --- #
    # ------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( ksteps  is None ): ksteps  = np.arange( config["arg_Iter"] )*config["arg_Step"]+config["arg_Init"]
    if ( datDir  is None ): datDir  = "{0}dat/".format( jobDir )
    if ( pngDir  is None ): pngDir  = "{0}png/".format( jobDir )
    if ( datFile is None ): datFile = OutDir + "EdotJAnalysis.dat"
    if ( pngFile is None ): pngFile = OutDir + "EdotJAnalysis.png"
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    if ( uorJ    is None ): uorJ    = "J"
    # ------------------------------- #
    # --- [2]  J.E を解析         --- #
    # ------------------------------- #
    if ( Flag__CallCalculater ):
        Nsteps = len( ksteps )
        if ( NewFile ): ( open( OutFile, "w" ) ).close()    # -- 新規ファイル -- #
        pgb.progressBar( iteration=0, total=Nsteps )        # -- 進捗表示     -- #
        for i,kstep in enumerate( ksteps ):
            ret = EdJ_calculater( job=job, jobDir=jobDir, kstep=kstep, config=config, uorJ=uorJ )
            pgb.progressBar( iteration=i+1, total=Nsteps )
            with open( OutFile, "ab" ) as f:
                np.savetxt( f, np.array( ret ).reshape( (1,-1) ) )
    # ------------------------------- #
    # --- [3] 返却 /  描画        --- #
    # ------------------------------- #
    #  -- [3-1] .png 書き出し     --  #
    if ( Flag__CallPlotter ):
        EdJ_plotter( datFile=OutFile, config=config, pngDir=pngDir, CnsFile=CnsFile )



# ======================================== #
# ===    J.Eを計算して表示             === #
# ======================================== #
def EdJ_calculator( job  =None, jobDir=None, kstep =None, \
                    uorJ =None, config=None ):
    # ------------------------ #
    # --- [1] 引数チェック --- #
    # ------------------------ #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"    .format( config["pic_jobDir"], job )
    if ( uorJ    is None ): uorJ    = "J"
    # -------------------------- #
    # --- [2] データ呼び出し --- #
    # -------------------------- #
    keys  = [ "Ex", "Ey", "Ez", "Jex", "Jey", "Jez", "Jix", "Jiy", "Jiz", "ne", "ni", "psi" ]
    config["cmp_nFilter"]    = 10
    config["cmp_LinearFilt"] = 0.2
    Data  = fpc.FetchPIC( job=job , jobDir=jobDir, kstep =kstep, keys=keys, config=config, FilterKeys=keys )
    xAxis = + Data["xAxis"]
    yAxis = + Data["yAxis"]
    flux  = + Data["psi"  ]
    if ( uorJ == "J" ):
        Fex   =   Data["Jex"]
        Fey   =   Data["Jey"]
        Fez   =   Data["Jez"]
        Fix   =   Data["Jix"]
        Fiy   =   Data["Jiy"]
        Fiz   =   Data["Jiz"]
    if ( uorJ == "u" ):
        Fex   = - Data["Jex"] * riv.robustInv( Data["ne"] )
        Fey   = - Data["Jey"] * riv.robustInv( Data["ne"] )
        Fez   = - Data["Jez"] * riv.robustInv( Data["ne"] )
        Fix   = + Data["Jix"] * riv.robustInv( Data["ni"] )
        Fiy   = + Data["Jiy"] * riv.robustInv( Data["ni"] )
        Fiz   = + Data["Jiz"] * riv.robustInv( Data["ni"] )
    # -------------------------- #
    # --- [3]  E.u , E.J     --- #
    # -------------------------- #
    #  -- [3-1] E.ue         --  #
    EdJepar  = 
    EdJeprp  = 
    EdJipar  = 
    EdJiprp  = 
    dotEp_e  = Fex*Data["Ex"] + Fez*Data["Ez"]
    dotEt_e  = Fey*Data["Ey"]
    dotEs_e  = dotEp_e + dotEt_e
    #  -- [3-2] E.ui         --  #
    dotEp_i  = Fix*Data["Ex"] + Fiz*Data["Ez"]
    dotEt_i  = Fiy*Data["Ey"]
    dotEs_i  = dotEp_i + dotEt_i
    # ------------------------- #
    # --- [4] 合計          --- #
    # ------------------------- #
    mask     = 1.0
    volm     = cvl.cylindricalVolume( zAxis=xAxis, rAxis=yAxis )
    sdotEp_e = np.sum( dotEp_e * volm * mask )
    sdotEt_e = np.sum( dotEt_e * volm * mask )
    sdotEs_e = np.sum( dotEs_e * volm * mask )
    sdotEp_i = np.sum( dotEp_e * volm * mask )
    sdotEt_i = np.sum( dotEt_e * volm * mask )
    sdotEs_i = np.sum( dotEs_e * volm * mask )
    return( [ fstep, ptime, sdotEp_e, sdotEt_e, sdotEs_e, sdotEp_i, sdotEt_i, sdotEs_i ] )


def timeEdJ_plotter( datFile=None )

    # ------------------------- #
    # --- [1] プロット      --- #
    # ------------------------- #
    cfs.configSettings( config=config, configType="plot1D_lateral" )
    cfs.configSettings( config=config, configType="FilterClear"    )
    FigName = pngFile.replace( "EdotJ", "EdotJe" )
    fig     = pl1.plot1D( FigName=FigName, config=config )
    fig.addPlot( xAxis= )
    fig.setLegend()
    fig.setAxis()
    fig.writeFigure()
    FigName = pngFile.replace( "EdotJ", "EdotJi" )
    print( "[time_EdJ] {0} is outputed...".format( FigName ) )


# ======================================== #
# ===           テスト実行             === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    OutDir  = "./png/t_EdJ/"
    for kstep in args["ksteps"]:
        time_EdJ( job  =args["job"], jobDir=args["jobDir"], \
                 kstep=kstep      , OutDir=OutDir          )
