import myStyle.cMap2D     as clm
import mPICsee.FetchHelix as fhx
import myStyle.LoadConfig as lcf
import strManipulate.keyInclude as stm

# --- pic から ヘリシティ 2次元カラーマップ描写 --- #
def hlx2cmp( job  =None, jobDir=None, kstep=None, OutFile=None, OutDir=None, \
             keys =None, config=None, Flag__FieldReturn=False ):
    # ------------------------- #
    # --- [1] 引数チェック  --- #
    # ------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"    .format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}png/"    .format( jobDir )
    if ( OutFile is None ): OutFile = "Helix2D.png"
    FieldVars = [ "uex", "uey", "uez", "uix", "uiy", "uiz", \
                  "wex", "wey", "wez", "wix", "wiy", "wiz", \
                  "Bx" , "By" , "Bz" , "Avx", "Avy", "Avz"  ]
    Flag__FieldReturn = stm.keyInclude( list1=keys, list2=FieldVars )
    OutFile           = OutDir + OutFile
    
    # ------------------------- #
    # --- [2] データ呼出    --- #
    # ------------------------- #
    Data  = fhx.FetchHelix( job   =job   , jobDir=jobDir, kstep =kstep, \
                            config=config, Flag__FieldReturn=Flag__FieldReturn )
    xAxis = Data["xAxis"]
    yAxis = Data["yAxis"]
    flux  = Data["flux" ]
    
    # ------------------------- #
    # --- [3]   プロット    --- #
    # ------------------------- #
    for key in keys:
        FigName = OutFile.replace( "Helix2D", "Helix2D_{0}_{1:08}".format( key, kstep )  )
        clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis, cMap=Data[key], Cntr=flux, \
                    FigName=FigName, config=config )

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
        p2c    = hlx2cmp( job  =job  , jobDir=jobDir, OutDir=OutDir, \
                          kstep=kstep, keys  =keys                   )
