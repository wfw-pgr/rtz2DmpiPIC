import myStyle.cMap2D         as clm
import mPICsee.FetchPIC       as fpc
import myStyle.LoadConfig     as lcf
import myStyle.configSettings as cfs

# --- pic から 2次元カラーマップ描写 --- #
def Epp2cmp( job  =None, jobDir=None, kstep=None, pngFile=None, pngDir=None, \
             config=None ):
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
    Fkeys =["psi", "Bx", "By", "Bz", "Ex", "Ey", "Ez" ]
    Ekeys =["vpara", "vprp1", "vprp2"]
    Data  = fpc.FetchPIC( job =job  , jobDir=jobDir, kstep=kstep, \
                          keys=Fkeys, config=config )
    xAxis = Data["xAxis"]
    yAxis = Data["yAxis"]
    import fLIB.fLIB__MFppVector as ppD
    Epp   = ppD.MFppVector( v1=Data["Ex"], v2=Data["Ey"], v3=Data["Ez"],\
                            r1=Data["Bx"], r2=Data["By"], r3=Data["Bz"] )
    # import fLIB.fLIB__MFppDecompose as ppD
    # unitB = ppD.MFppDecompose( r1=Data["Bx"], r2=Data["By"], r3=Data["Bz"] )
    # Epp = {}
    # Epp["vpara"] = unitB["vparax"]*Data["Ex"] + unitB["vparay"]*Data["Ey"] + unitB["vparaz"]*Data["Ez"]
    # Epp["vprp1"] = unitB["vprp1x"]*Data["Ex"] + unitB["vprp1y"]*Data["Ey"] + unitB["vprp1z"]*Data["Ez"]
    # Epp["vprp2"] = unitB["vprp2x"]*Data["Ex"] + unitB["vprp2y"]*Data["Ey"] + unitB["vprp2z"]*Data["Ez"]
    
    # ------------------------------ #
    # --- [3] プロット           --- #
    # ------------------------------ #
    cfs.configSettings( config=config, configType="cMap2D_thesis" )
    config["cmp_ColorTable"] = "bwr"
    config["cmp_AutoLevel"]  = True
    config["cmp_MaxMin"]     = [-0.2,+0.2]
    for key in Ekeys:
        FigName = pngFile.replace( "Field2D", "Field2D_{0}_{1:08}".format( key.replace("v","E"), kstep )  )
        clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis, cMap=Epp[key], Cntr=Data["psi"], \
                    FigName=FigName, config=config )

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    pngDir = "./png/"
    for kstep in args["ksteps"]:
        print( " [Epp2cmp] ( job , kstep ) = ( {0} , {1} )".format( args["job"], kstep ) )
        p2c    = Epp2cmp( job  =args["job"], jobDir=args["jobDir"], pngDir=pngDir, \
                          kstep=kstep                       )
