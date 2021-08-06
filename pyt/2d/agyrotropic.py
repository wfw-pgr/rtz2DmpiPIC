import myStyle.cMap2D           as clm
import mPICsee.FetchPIC         as fpc
import myStyle.LoadConfig       as lcf
import strManipulate.keyInclude as stm
import myStyle.configSettings   as cfs
import fLIB.fLIB__robustInv2D   as rI2

# --- pic から 2次元カラーマップ描写 --- #
def agyrotropic( job  =None, jobDir=None, kstep=None, OutFile=None, OutDir=None, \
                 keys =None, config=None ):
    # ------------------------------ #
    # --- [1] 引数チェック       --- #
    # ------------------------------ #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"    .format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}png/"    .format( jobDir )
    if ( OutFile is None ): OutFile = "Temperature.png".format( kstep  )
    OutFile = OutDir + OutFile
    
    # ------------------------------ #
    # --- [2] データ読込         --- #
    # ------------------------------ #
    Pkeys = [ "pen1", "pen2", "pet1", "psi" ]
    config["pic_Temperature"] = False
    Data        = fpc.FetchPIC( job =job  , jobDir=jobDir, kstep =kstep, \
                                keys=Pkeys, config=config )
    xAxis       = + Data["xAxis"]
    yAxis       = + Data["yAxis"]
    flux        = + Data["psi"  ]
    pPerp       = 0.5*( Data["pen1"] + Data["pen2"] )
    pPerpInv    = rI2.robustInv2D( Data=pPerp  )
    pTotal      = ( Data["pen1"]+Data["pen2"]+Data["pet1"] ) / 3.0
    pTotalInv   = rI2.robustInv2D( Data=pTotal )
    agyro1      = ( Data["pen1"]-Data["pen2"] ) * pPerpInv
    agyro2      = ( pPerp-Data["pet1"]        ) * pTotalInv
    
    # ------------------------------ #
    # --- [3] プロット           --- #
    # ------------------------------ #
    cfs.configSettings( config=config, configType="cMap2D_def" )
    config["cmp_ColorTable"] = "hot"
    config["cmp_AutoLevel"]  = True
    config["cmp_MaxMin"]     = [0.0,0.12]
    config["MinimalOut"]     = True
    FigName = OutFile.replace( ".png", "_agyro1_{0:08}.png".format( kstep )  )
    clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis, cMap=agyro1, Cntr=flux, \
                FigName=FigName, config=config )
    FigName = OutFile.replace( ".png", "_agyro2_{0:08}.png".format( kstep )  )
    clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis, cMap=agyro2, Cntr=flux, \
                FigName=FigName, config=config )

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    OutDir = "tmp/"
    for kstep in args["ksteps"]:
        print( " [Tei2cmp] ( job , kstep ) = ( {0} , {1} )".format( args["job"], kstep ), end=""  )
        p2c    = agyrotropic( job  =args["job"], jobDir=args["jobDir"], OutDir=OutDir, \
                          kstep=kstep      , keys  =args["key"]                    )
        print( "\t [Done]" )
