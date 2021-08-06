import myStyle.cMap2D         as clm
import myStyle.LoadConfig     as lcf
import myStyle.configSettings as cfs
import mPICsee.driftField     as dfd
import mPICsee.FetchPIC       as fpc

# ======================================== #
# ===  ドリフト流体場を表示            === #
# ======================================== #
def driftF2cmp( job    =None, jobDir=None, kstep=None, OutFile=None, OutDir=None, \
                RawData=None, config=None ):
    # ------------------------------ #
    # --- [1] 引数チェック       --- #
    # ------------------------------ #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"       .format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}png/"       .format( jobDir )
    if ( OutFile is None ): OutFile = "{0}Field2D.png".format( OutDir, kstep )
    
    # ------------------------------ #
    # --- [2] データ読込         --- #
    # ------------------------------ #
    Fkeys = ["psi","Bx"]
    dkeys = ["vExBx", "vExBy", "vExBz", "vgrBx", "vgrBy", "vgrBz", \
             "vcvtx", "vcvty", "vcvtz", "vmagx", "vmagy", "vmagz", "vplrx", "vplry", "vplrz"  ]
    Data  = fpc.FetchPIC( job =job  , jobDir=jobDir, kstep =kstep, \
                          keys=Fkeys, config=config, RawData=RawData )
    drift = dfd.driftField( job=job, jobDir=jobDir, kstep=kstep, config=config, RawData=RawData )
    xAxis = Data["xAxis"]
    yAxis = Data["yAxis"]
    flux  = Data["psi"  ]
    
    # ------------------------------ #
    # --- [3] プロット           --- #
    # ------------------------------ #
    cfs.configSettings( config=config, configType="cMap2D_thesis" )
    config["cmp_AutoLevel"]  = False
    config["cmp_MaxMin"]     = [-0.1,0.1]
    config["cmp_ColorTable"] = "jet"
    for key in dkeys:
        cMap    = drift[ key ]
        FigName = OutFile.replace( ".png", "_{0}_{1:08}.png".format( key, kstep )  )
        clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis, cMap=cMap, Cntr=flux, \
                    FigName=FigName, config=config )


# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    OutDir = "./png/driftF2cmp/"
    for kstep in args["ksteps"]:
        print(" [driftF2cmp] ( job, kstep ) = ( {0} , {1} )".format( args["job"], kstep ) )
        p2c = driftF2cmp( job=args["job"], jobDir=args["jobDir"], OutDir=OutDir, kstep=kstep )
