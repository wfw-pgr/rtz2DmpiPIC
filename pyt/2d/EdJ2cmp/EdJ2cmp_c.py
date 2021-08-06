import myStyle.cMap2D            as clm
import mPICsee.FetchPIC          as fpc
import myStyle.LoadConfig        as lcf
import myStyle.configSettings    as cfs
import fLIB.fLIB__RefppDecompose as ppD
import mPICsee.EstidDecompose    as est

# ======================================== #
# ===    J.Eを計算して表示             === #
# ======================================== #
def EdJ2cmp( job  =None, jobDir=None, kstep  =None, OutFile=None, OutDir=None, \
             keys =None, config=None, RawData=False ):
    # ------------------------ #
    # --- [1] 引数チェック --- #
    # ------------------------ #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    jobDir = "../" + jobDir
    if ( OutDir  is None ): OutDir  = "{0}png/".format( jobDir )
    if ( OutFile is None ): OutFile = "{0}EdJ2cmp__.png".format( OutDir  )
    # -------------------------- #
    # --- [2] データ呼び出し --- #
    # -------------------------- #
    fkeys  = [ "Bx", "By", "Bz", "Ex", "Ey", "Ez", "Jex", "Jey", "Jez", "Jix", "Jiy", "Jiz", "psi" ]
    config["cmp_nFilter"]    = 20
    config["cmp_LinearFilt"] = 0.3
    Data  = fpc.FetchPIC( job=job , jobDir=jobDir, kstep =kstep, keys=fkeys, config=config, FilterKeys=fkeys, RawData=RawData )
    xAxis = + Data["xAxis"]
    yAxis = + Data["yAxis"]
    flux  = + Data["psi"  ]
    EppD  = ppD.RefppDecompose( v1=Data["Ex"], v2=Data["Ey"], v3=Data["Ez"], r1=Data["Bx"], r2=Data["By"], r3=Data["Bz"] )
    EsiD  = est.EstidDecompose( job=job, jobDir=jobDir, kstep=kstep, RawData=RawData )
    # -------------------------- #
    # --- [3]  E.u , E.J     --- #
    # -------------------------- #
    #  -- [3-1] E.ue         --  !
    Data["pol_e"] = Data["Jex"]*Data["Ex"] + Data["Jez"]*Data["Ez"]
    Data["tor_e"] = Data["Jey"]*Data["Ey"]
    Data["prl_e"] = Data["Jex"]*EppD["vparax"] + Data["Jey"]*EppD["vparay"] + Data["Jez"]*EppD["vparaz"]
    Data["prp_e"] = Data["Jex"]*EppD["vprp1x"] + Data["Jey"]*EppD["vprp1y"] + Data["Jez"]*EppD["vprp1z"]
    Data["stE_e"] = Data["Jex"]*EsiD["Estx"  ] + Data["Jey"]*EsiD["Esty"  ] + Data["Jez"]*EsiD["Estz"  ]
    Data["idE_e"] = Data["Jex"]*EsiD["Eidx"  ] + Data["Jey"]*EsiD["Eidy"  ] + Data["Jez"]*EsiD["Eidz"  ]
    Data["EdJ_e"] = Data["pol_e"] + Data["tor_e"]
    #  -- [3-2] E.ui         --  !
    Data["pol_i"] = Data["Jix"]*Data["Ex"] + Data["Jiz"]*Data["Ez"]
    Data["tor_i"] = Data["Jiy"]*Data["Ey"]
    Data["prl_i"] = Data["Jix"]*EppD["vparax"] + Data["Jiy"]*EppD["vparay"] + Data["Jiz"]*EppD["vparaz"]
    Data["prp_i"] = Data["Jix"]*EppD["vprp1x"] + Data["Jiy"]*EppD["vprp1y"] + Data["Jiz"]*EppD["vprp1z"]
    Data["stE_i"] = Data["Jix"]*EsiD["Estx"  ] + Data["Jiy"]*EsiD["Esty"  ] + Data["Jiz"]*EsiD["Estz"  ]
    Data["idE_i"] = Data["Jix"]*EsiD["Eidx"  ] + Data["Jiy"]*EsiD["Eidy"  ] + Data["Jiz"]*EsiD["Eidz"  ]
    Data["EdJ_i"] = Data["pol_i"] + Data["tor_i"]
    # ------------------------- #
    # --- [4] プロット      --- #
    # ------------------------- #
    cfs.configSettings( config=config, configType="cMap2D_thesis" )
    cfs.configSettings( config=config, configType="NoAxis" )    
    config["cmp_ColorTable"] = "jet"
    config["cmp_AutoLevel"]  = False
    config["cmp_MaxMin"]     = [-0.01,0.01]
    config["clb_sw"]         = False
    config["cmp_yAutoRange"] = False
    config["cmp_yRange"]     = [-20.0,20.0]
    pkeys = [ "pol_e", "tor_e", "prp_e", "prl_e", "stE_e", "idE_e", "EdJ_e", \
              "pol_i", "tor_i", "prp_i", "prl_i", "stE_i", "idE_i", "EdJ_i"  ]
    pkeys = [ "pol_i"]
    for key in pkeys:
        FigName = OutFile.replace( "__.png", "_{0}_{1:08}.png".format( key, kstep ) )
        clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis , cMap=Data[key], Cntr=flux, \
                    FigName=FigName, config=config )


# ======================================== #
# ===           テスト実行             === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    OutDir  = "./png/"
    for kstep in args["ksteps"]:
        EdJ2cmp( job  =args["job"], jobDir=args["jobDir"], \
                 kstep=kstep      , OutDir=OutDir          )
