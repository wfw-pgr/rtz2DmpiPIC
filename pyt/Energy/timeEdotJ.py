import numpy                         as np
import mPICsee.FetchPIC              as fpc
import mPICsee.TimeUnitConvert       as tuc
import myUtils.progressBar           as pgb
import myStyle.plot1D                as pl1
import myStyle.LoadConfig            as lcf
import myStyle.configSettings        as cfs
import myBasicAlgs.robustInv         as riv
import myBasicAlgs.cylindricalVolume as cvl
import myAnalysis.createMask         as msk
import myAnalysis.outletMask         as olt
import myConvert.arr2dct             as a2d

# ======================================== #
# ===    J.Eを計算して表示             === #
# ======================================== #
def EdJ2cmp( job    =None , jobDir =None, ksteps=None , OutFile=None, OutDir=None, pngDir=None, \
             uorJ   =None , retSize=12  , config=None , CnsFile=None, \
             NewFile=True , Flag__CallCalculater=True , Flag__CallPlotter=True, Flag__timeIntegration=True ):
    # ------------------------------- #
    # --- [1] 引数チェック        --- #
    # ------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( ksteps  is None ): ksteps  = np.arange( config["arg_Iter"] )*config["arg_Step"]+config["arg_Init"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}dat/".format( jobDir )
    if ( pngDir  is None ): pngDir  = "{0}png/".format( jobDir )
    if ( OutFile is None ): OutFile = "{0}t-EdotJ.dat".format( OutDir )
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    if ( uorJ    is None ): uorJ    = "J"
    # ------------------------------- #
    # --- [2] E.J / E.u 計算      --- #
    # ------------------------------- #
    if ( Flag__CallCalculater ):
        Nsteps = len( ksteps )
        if ( NewFile ): ( open( OutFile, "w" ) ).close()    # -- 新規ファイル -- #
        pgb.progressBar( iteration=0, total=Nsteps )        # -- 進捗表示     -- #
        for i,kstep in enumerate( ksteps ):
            ret = EdJ_calculator( job =job , jobDir=jobDir, kstep=kstep, \
                                  uorJ=uorJ, config=config )
            pgb.progressBar( iteration=i+1, total=Nsteps )
            with open( OutFile, "ab" ) as f:
                np.savetxt( f, np.array( ret ).reshape( (1,retSize) ) )
    # ------------------------------- #
    # --- [3] 返却 /  描画        --- #
    # ------------------------------- #
    #  -- [3-1] .png 書き出し     --  #
    if ( Flag__CallPlotter ):
        EdJ_plotter( datFile=OutFile, config=config, pngDir=pngDir, CnsFile=CnsFile, \
                     Flag__timeIntegration=Flag__timeIntegration )


def EdJ_calculator( job=None, jobDir=None, kstep=None, uorJ=None, config=None ):
    # -------------------------- #
    # --- [1] データ呼び出し --- #
    # -------------------------- #
    keys  = [ "Ex", "Ey", "Ez", "Jex", "Jey", "Jez", "Jix", "Jiy", "Jiz", "ne", "ni", "psi" ]
    config["phys_TimeUnit"] = "wce"
    Data  = fpc.FetchPIC( job=job , jobDir=jobDir, kstep =kstep, keys=keys, config=config )
    if ( uorJ == "J" ):
        Fex, Fey, Fez = Data["Jex"], Data["Jey"], Data["Jez"]
        Fix, Fiy, Fiz = Data["Jix"], Data["Jiy"], Data["Jiz"]
    if ( uorJ == "u" ):
        neInv         = riv.robustInv( Data["ne"] )
        niInv         = riv.robustInv( Data["ni"] )
        Fex, Fey, Fez = - Data["Jex"]*neInv, - Data["Jey"]*neInv, - Data["Jez"]*neInv
        Fix, Fiy, Fiz = - Data["Jix"]*niInv, - Data["Jiy"]*niInv, - Data["Jiz"]*niInv
    # -------------------------- #
    # --- [2]  E.u , E.J     --- #
    # -------------------------- #
    #  -- [2-1] E.ue         --  !
    fstep     = float( kstep )
    ptime     = Data["time"]
    volume    = cvl.cylindricalVolume( zAxis=Data["xAxis"], rAxis=Data["yAxis"] )
    masks     = msk.createMask( Flux=Data["psi"], Flag__privateFluxMode=True )
    omask     = olt.outletMask( Flux=Data["psi"], psiSurpass=0.1 )
    dotE_e    = Fex*Data["Ex"] + Fey*Data["Ey"] + Fez*Data["Ez"]
    dotE_i    = Fix*Data["Ex"] + Fiy*Data["Ey"] + Fiz*Data["Ez"]
    dotEe_f   = np.sum( dotE_e * volume * masks["fMask3" ] )
    dotEe_o1  = np.sum( dotE_e * volume * omask["outlet1"] )
    dotEe_o2  = np.sum( dotE_e * volume * omask["outlet2"] )
    dotEe_w   = np.sum( dotE_e * volume )
    dotEe_d   = dotEe_w - dotEe_f
    dotEi_f   = np.sum( dotE_i * volume * masks["fMask3" ] )
    dotEi_o1  = np.sum( dotE_e * volume * omask["outlet1"] )
    dotEi_o2  = np.sum( dotE_e * volume * omask["outlet2"] )
    dotEi_w   = np.sum( dotE_i * volume )
    dotEi_d   = dotEi_w - dotEi_f
    return( [ fstep, ptime, \
              dotEe_f , dotEi_f , dotEe_d , dotEi_d , dotEe_w, dotEi_w, \
              dotEe_o1, dotEi_o1, dotEe_o2, dotEi_o2 ] )
    

# --- .png 書き出し --- #
def EdJ_plotter( datFile=None, config=None, pngDir=None, CnsFile=None, Flag__timeIntegration=True ):
    # ---------------------------------------- #
    # --- [1] プロット準備                  -- #
    # ---------------------------------------- #
    #  -- [1-1] データ読込                 --  #
    Data  = np.transpose( np.loadtxt( datFile ) )
    keys  = ["fstep","ptime","dotEe_f","dotEi_f","dotEe_d","dotEi_d","dotEe_w","dotEi_w", \
             "dotEe_o1","dotEi_o1","dotEe_o2","dotEi_o2"]
    Data  = a2d.arr2dct( Data=Data, keys=keys )
    ptime = tuc.TimeUnitConvert( Data=Data["ptime"], unit="wpi", CnsFile=CnsFile )
    if ( Flag__timeIntegration ):
        ikeys = ["dotEe_f","dotEi_f","dotEe_d","dotEi_d","dotEe_w","dotEi_w","dotEe_o1","dotEi_o1","dotEe_o2","dotEi_o2"]
        import myBasicAlgs.integrate1D as it1
        for key in ikeys:
            Data[key] = it1.integrate1D( xAxis=Data["ptime"], yAxis=Data[key], initial=0.0 )
            
    # ---------------------------------------- #
    # --- [2] プロット                     --- #
    # ---------------------------------------- #
    #  -- [1-1] コンフィグ                 --  #
    cfs.configSettings( config=config, configType="plot1D_lateral" )
    cfs.configSettings( config=config, configType="FilterClear"    )
    config["plt_xAutoRange"]  = False
    config["plt_xRange"]      = [0.,120.]
    config["plt_yAutoRange"]  = True
    config["yMajor_Nticks"]   = 7
    config["plt_GaussFilt"]   = 0.0

    #  -- [1-2] Electron ( 電子 ) -- #
    OutFile = pngDir + "t-EdotJe.png"
    fig = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=ptime, yAxis=Data["dotEe_f" ] )
    fig.addPlot( xAxis=ptime, yAxis=Data["dotEe_o1"] )
    fig.addPlot( xAxis=ptime, yAxis=Data["dotEe_o2"] )
    fig.addPlot( xAxis=ptime, yAxis=Data["dotEe_d" ] )
    fig.addPlot( xAxis=ptime, yAxis=Data["dotEe_w" ] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    #  -- [1-3]  ion   ( イオン ) -- #
    OutFile = pngDir + "t-EdotJi.png"
    fig = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=ptime, yAxis=Data["dotEi_f" ] )
    fig.addPlot( xAxis=ptime, yAxis=Data["dotEi_o1"] )
    fig.addPlot( xAxis=ptime, yAxis=Data["dotEi_o2"] )
    fig.addPlot( xAxis=ptime, yAxis=Data["dotEi_d" ] )
    fig.addPlot( xAxis=ptime, yAxis=Data["dotEi_w" ] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()

# ======================================== #
# ===           テスト実行             === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    pngDir  = "./timeEdotJ/"
    EdJ2cmp( job   =args["job"]   , jobDir=args["jobDir"], \
             ksteps=args["ksteps"], pngDir=pngDir          )
