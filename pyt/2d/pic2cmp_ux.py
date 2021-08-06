import numpy as np
import myStyle.cMap2D         as clm
import mPICsee.FetchPIC       as fpc
import myStyle.LoadConfig     as lcf
import myStyle.configSettings as cfs
import mPICsee.SpeedUnitConvert  as suc
import mPICsee.LengthUnitConvert as luc

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
    Dkeys       = ["psi", "Bx", "By", "Bz", "uix", "uiy" ]
    Data        = fpc.FetchPIC( job =job  , jobDir=jobDir, kstep =kstep, \
                                keys=Dkeys, config=config )
    xAxis       = + Data["xAxis"]
    yAxis       = + Data["yAxis"]
    flux        = + Data["psi"  ]
    absB        = np.mean( Data["Bx"]**2 + Data["By"]**2 + Data["Bz"]**2 )
    import myBasicAlgs.robustInv as inv
    Bcoef       = np.sqrt( inv.robustInv( absB ) )
    # ------------------------------ #
    # --- [3] プロット           --- #
    # ------------------------------ #
    cfs.configSettings( config=config, configType="cMap2D_thesis" )
    config["cmp_ColorTable"] = "bwr"
    config["FigSize"]        = (6,4)
    config["cmp_AutoLevel"]  = False
    config["cmp_MaxMin"]     = [-0.8,+0.8]
    config["cmp_xAutoRange"] = False
    config["cmp_xRange"]     = [0.0,+40.0]
    config["cmp_yAutoRange"] = False
    config["cmp_yRange"]     = [-12.0,+12.0]
    config["xMajor_FontSize"]= 16
    config["yMajor_FontSize"]= 16
    config["clb_Nlabel"]     = 2
    config["clb_FontSize"]   = 16
    config["clb_Orientation"]= "vertical"
    config["cmp_Position"]   = [0.15,0.15,0.85,0.85]
    config["clb_Position"]   = [0.88,0.185,0.92,0.815]
    config["phys_SpeedUnit"] = "vAi"
    ret         = luc.LengthUnitConvert( xAxis=xAxis, yAxis=yAxis, job=job, jobDir=jobDir, config=config )
    xAxis,yAxis = ret["xAxis"], ret["yAxis"]
    Data["uix"] = suc.SpeedUnitConvert( Data=Data["uix"], job=job, jobDir=jobDir, config=config ) * Bcoef
    Data["uiy"] = suc.SpeedUnitConvert( Data=Data["uiy"], job=job, jobDir=jobDir, config=config ) * Bcoef
    ukeys       = ["uix","uiy"]
    for key in ukeys:
        FigName = OutFile.replace( ".png", "_{0}_{1:08}.png".format( key, kstep )  )
        clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis, cMap=Data[ key ], Cntr=flux, \
                    FigName=FigName, config=config )

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    job    = args["job"]
    jobDir = args["jobDir"]
    kstep  = args["kstep"]
    OutDir = "./png/uField2D/"
    for kstep in args["ksteps"]:
        print( " [pic2cmp] ( job , kstep ) = ( {0} , {1:8} )".format( job, kstep ), end=""  )
        p2c    = pic2cmp( job  =job  , jobDir=jobDir, OutDir=OutDir, \
                          kstep=kstep )
        print("\t [Done]")
