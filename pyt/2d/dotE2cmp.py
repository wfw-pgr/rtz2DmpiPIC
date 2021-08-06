import myStyle.cMap2D        as clm
import mPICsee.FetchPIC      as fpc
import myStyle.LoadConfig    as lcf
import myBasicAlgs.robustInv as riv
import myUtils.LoadConst     as lcn
import myStyle.configSettings as cfs

# ======================================== #
# ===    J.Eを計算して表示             === #
# ======================================== #
def dotE2cmp( job  =None, jobDir=None, kstep=None, OutFile=None, OutDir=None, CnsFile=None, \
              keys =None, config=None ):
    # ------------------------ #
    # --- [1] 引数チェック --- #
    # ------------------------ #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"             .format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}png/"             .format( jobDir )
    if ( OutFile is None ): OutFile = "dotE2D.png"          .format( kstep  )
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    OutFile = OutDir + OutFile

    # -------------------------- #
    # --- [2] データ呼び出し --- #
    # -------------------------- #
    const  = lcn.LoadConst( InpFile=CnsFile )
    keys  = [ "Ex", "Ey", "Ez", "Jex", "Jey", "Jez", "Jix", "Jiy", "Jiz", "ne", "ni", "psi" ]
    Data  = fpc.FetchPIC( job =job , jobDir=jobDir, kstep  =kstep, \
                          keys=keys, config=config, RawData=False  )
    xAxis = + Data["xAxis"]
    yAxis = + Data["yAxis"]
    flux  = + Data["psi"  ]
    Data["uex"] = - Data["Jex"] * riv.robustInv( Data["ne"] )
    Data["uey"] = - Data["Jey"] * riv.robustInv( Data["ne"] )
    Data["uez"] = - Data["Jez"] * riv.robustInv( Data["ne"] )
    Data["uix"] = + Data["Jix"] * riv.robustInv( Data["ni"] )
    Data["uiy"] = + Data["Jiy"] * riv.robustInv( Data["ni"] )
    Data["uiz"] = + Data["Jiz"] * riv.robustInv( Data["ni"] )

    # -------------------------- #
    # --- [3]  E.u           --- #
    # -------------------------- #
    #  -- [3-1] E.ue         --  !
    uedotEp = Data["uex"]*Data["Ex"] + Data["uez"]*Data["Ez"]
    uedotEt = Data["uey"]*Data["Ey"]
    uedotE  = uedotEp + uedotEt
    #  -- [3-2] E.ui         --  !
    uidotEp = Data["uix"]*Data["Ex"] + Data["uiz"]*Data["Ez"]
    uidotEt = Data["uiy"]*Data["Ey"]
    uidotE  = uidotEp + uidotEt
    
    # -------------------------- #
    # --- [4]  E.J           --- #
    # -------------------------- #
    #  -- [4-1] E.Je         --  !
    JedotEp = Data["Jex"]*Data["Ex"] + Data["Jez"]*Data["Ez"]
    JedotEt = Data["Jey"]*Data["Ey"]
    JedotE  = JedotEp + JedotEt
    #  -- [4-2] E.Ji         --  !
    JidotEp = Data["Jix"]*Data["Ex"] + Data["Jiz"]*Data["Ez"]
    JidotEt = Data["Jiy"]*Data["Ey"]
    JidotE  = JidotEp + JidotEt
    
    # ------------------------- #
    # --- [5] プロット      --- #
    # ------------------------- #
    cfs.configSettings( config=config, configType="cMap2D_thesis" )
    config["cmp_AutoLevel"] = True
    cMaps  = { "uedotEp":uedotEp, "uedotEt":uedotEt, "uedotE":uedotE, "uidotEp":uidotEp, "uidotEt":uidotEt, "uidotE":uidotE,\
               "JedotEp":JedotEp, "JedotEt":JedotEt, "JedotE":JedotE, "JidotEp":JidotEp, "JidotEt":JidotEt, "JidotE":JidotE }
    pKeys  = cMaps.keys()
    for key in pKeys:
        FigName = OutFile.replace( ".png", "_{0}_{1:08}.png".format( key, kstep )  )
        clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis , cMap=cMaps[key], Cntr=flux, \
                    FigName=FigName, config=config )
        print( "[dotE2cmp] {0} is outputed...".format( FigName ) )


# ------------------------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    OutDir  = "./dotE2cmp/"
    for kstep in args["ksteps"]:
        dotE2cmp( job  =args["job"], jobDir=args["jobDir"], \
                  kstep=kstep      , OutDir=OutDir          )
