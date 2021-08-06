import myStyle.cMap2D         as clm
import mPICsee.FetchPIC       as fpc
import myStyle.LoadConfig     as lcf
import myStyle.configSettings as cfs

# ======================================== #
# === pic から 2次元カラーマップ描写   === #
# ======================================== #
def pic2cmp( job  =None, jobDir=None, kstep=None, pngFile=None, pngDir=None, \
             keys =None, config=None ):
    # ------------------------------ #
    # --- [1] 引数チェック       --- #
    # ------------------------------ #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"       .format( config["pic_jobDir"], job )
    if ( pngDir  is None ): pngDir  = "{0}png/"       .format( jobDir )
    if ( pngFile is None ): pngFile = "{0}Field2D.png".format( pngDir  )
    # ------------------------------ #
    # --- [2] データ読込         --- #
    # ------------------------------ #
    if ( not( "psi" in keys )): Fkeys = keys + ["psi"]
    Data  = fpc.FetchPIC( job =job  , jobDir=jobDir, kstep=kstep, \
                          keys=Fkeys, config=config )
    xAxis = Data["xAxis"]
    yAxis = Data["yAxis"]
    # ------------------------------ #
    # --- [3] プロット           --- #
    # ------------------------------ #
    cfs.configSettings( config=config, configType="cMap2D_thesis" )
    config["cmp_ColorTable"] = "bwr"
    config["cmp_AutoLevel"]  = True
    config["cmp_MaxMin"]     = [-0.5,+0.5]
    for key in keys:
        FigName = pngFile.replace( ".png", "_{0}_{1:08}.png".format( key, kstep )  )
        clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis, cMap=Data[key], Cntr=Data["psi"], \
                    FigName=FigName, config=config )

        
# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    pngDir = "./png/"
    for kstep in args["ksteps"]:
        p2c    = pic2cmp( job  =args["job"], jobDir=args["jobDir"], pngDir=pngDir, \
                          kstep=kstep      , keys  =args["key"]                    )
