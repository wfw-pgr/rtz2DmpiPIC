import myStyle.cMap2D         as clm
import mPICsee.FetchPIC       as fpc
import myStyle.LoadConfig     as lcf
import myBasicAlgs.robustInv  as riv
import myStyle.configSettings as cfs

# ======================================== #
# ===    J.Eを計算して表示             === #
# ======================================== #
def EdJ2cmp( job  =None, jobDir=None, kstep=None, OutFile=None, OutDir=None, \
             uorJ =None, keys =None, config=None ):
    # ------------------------ #
    # --- [1] 引数チェック --- #
    # ------------------------ #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"    .format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}png/"    .format( jobDir )
    if ( OutFile is None ): OutFile = "t-EdotJ.png".format( kstep  )
    if ( uorJ    is None ): uorJ    = "J"
    OutFile = OutDir + OutFile
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
    #  -- [3-1] E.ue         --  !
    dotEp_e = Fex*Data["Ex"] + Fez*Data["Ez"]
    dotEt_e = Fey*Data["Ey"]
    dotEs_e = dotEp_e + dotEt_e
    #  -- [3-2] E.ui         --  !
    dotEp_i = Fix*Data["Ex"] + Fiz*Data["Ez"]
    dotEt_i = Fiy*Data["Ey"]
    dotEs_i = dotEp_i + dotEt_i
    
    # ------------------------- #
    # --- [5] プロット      --- #
    # ------------------------- #
    cfs.configSettings( config=config, configType="cMap2D_thesis" )
    config["cmp_ColorTable"] = "jet"
    config["cmp_AutoLevel"]  = True
    cMaps  = { "dotEp_e":dotEp_e, "dotEt_e":dotEt_e, "dotEs_e":dotEs_e, "dotEp_i":dotEp_i, "dotEt_i":dotEt_i, "dotEs_i":dotEs_i }
    pKeys  = cMaps.keys()
    for key in pKeys:
        FigName = OutFile.replace( ".png", "_{0}_{1:08}.png".format( key, kstep )  )
        clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis , cMap=cMaps[key], Cntr=flux, \
                    FigName=FigName, config=config )
        print( "[EdJ2cmp] {0} is outputed...".format( FigName ) )

# ======================================== #
# ===           テスト実行             === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    OutDir  = "./EdJ2cmp/"
    for kstep in args["ksteps"]:
        EdJ2cmp( job  =args["job"], jobDir=args["jobDir"], \
                 kstep=kstep      , OutDir=OutDir          )
