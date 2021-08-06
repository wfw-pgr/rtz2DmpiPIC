import myStyle.cMap2D         as clm
import mPICsee.FetchPIC       as fpc
import myStyle.LoadConfig     as lcf
import myStyle.configSettings as cfs
import Filter.LinearFilter2D  as flt

# --- pic から 2次元カラーマップ描写 --- #
def pic2cmp( job  =None, jobDir=None, kstep=None, OutFile=None, OutDir=None, \
             keys =None, config=None ):
    # ------------------------------ #
    # --- [1] 引数チェック       --- #
    # ------------------------------ #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"    .format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}png/"    .format( jobDir )
    if ( OutFile is None ): OutFile = "Field2D.png".format( kstep  )
    OutFile = OutDir + OutFile
    
    # ------------------------------ #
    # --- [2] データ読込         --- #
    # ------------------------------ #
    cmkey  = "Jy"
    uxkey  = "uix"
    uykey  = "uiz"
    Fkeys  = ["psi",cmkey, uxkey, uykey ]
    Data   = fpc.FetchPIC( job =job  , jobDir=jobDir, kstep =kstep, \
                           keys=Fkeys, config=config )
    xAxis  = Data["xAxis"]
    yAxis  = Data["yAxis"]
    flux   = Data["psi"]
    cMap   = Data[cmkey]
    xvec   = Data[uxkey]
    yvec   = Data[uykey]
    
    # ------------------------------ #
    # --- [3] プロット           --- #
    # ------------------------------ #
    cfs.configSettings( config=config, configType="cMap2D_thesis" )
    cfs.configSettings( config=config, configType="NoAxis"        )
    config["clb_sw"]         = False
    config["cmp_ColorTable"] = "bwr"
    config["cmp_AutoLevel"]  = False
    config["cmp_MaxMin"]     = [-0.2,+0.2]
    config["clb_Nlabel"]     = 0
    config["cmp_LinearFilt"] = 0.1
    cMap = flt.LinearFilter2D( Data=cMap, alpha=config["cmp_LinearFilt"], \
                               x1Axis=xAxis, x2Axis=yAxis, coordinate="rtz" )

    FigName = OutFile.replace( ".png", "_{0}_{1:08}.png".format( cmkey, kstep )  )
    clm.cMap2D( xAxis=yAxis, yAxis=xAxis, cMap=cMap      , Cntr  =flux, \
                xvec =xvec , yvec =yvec , FigName=FigName, config=config )
    

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    job    = args["job"]
    jobDir = args["jobDir"]
    kstep  = args["kstep"]
    keys   = args["key"]
    OutDir = "./png/"
    for kstep in args["ksteps"]:
        print( " [pic2cmp] ( job , kstep ) = ( {0} , {1:8} )".format( job, kstep ), end=""  )
        p2c    = pic2cmp( job  =job  , jobDir=jobDir, OutDir=OutDir, \
                          kstep=kstep, keys  =keys                   )
        print("\t [Done]")
