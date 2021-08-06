import numpy              as np
import myStyle.cMap2D     as clm
import mPICsee.FetchPIC   as fpc
import myStyle.LoadConfig as lcf
import myUtils.LoadConst  as lcn
import fLIB.fLIB__rct2rct as r2r
import myStyle.configSettings as cfs

# --- pic から 2次元カラーマップ描写 --- #
def grd2cmp( job    =None, jobDir=None, kstep=None, OutFile=None, OutDir=None, \
             CnsFile=None, dmnFile=None, config=None, Flag__PELine=True ):
    # ------------------------------ #
    # --- [1] 引数チェック       --- #
    # ------------------------------ #
    if ( config  is None ): config   = lcf.LoadConfig()
    if ( job     is None ): job      = config["pic_job"]
    if ( kstep   is None ): kstep    = config["pic_kstep"]
    if ( jobDir  is None ): jobDir   = "{0}{1}/"            .format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir   = "{0}png/"            .format( jobDir )
    if ( dmnFile is None ): dmnFile  = "{0}bin/kstep{1:08}/ijDomain{1:08}.dat".format( jobDir, kstep )
    if ( OutFile is None ): OutFile  = "Division_{0:08}.png".format( kstep  )
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    const    = lcn.LoadConst( InpFile=CnsFile )
    OutFile  = OutDir + OutFile
    ijDomain = np.array( np.loadtxt( dmnFile ), dtype=np.int32 )
    xAxiso   = np.linspace( -0.5*const["x1Leng"], +0.5*const["x1Leng"],     const["LI"]  )
    yAxiso   = np.linspace( const["x2Min"], const["x2Min"]+const["x2Leng"], const["LJ"]  )
    PEnum    = np.zeros( (const["LJ"], const["LI"]) )
    for i in range(16):
        PEnum[:,( ijDomain[i,5] ):( ijDomain[i,6]+1 )] = ijDomain[i,0]
    if ( config["pic_DownConvert"] ):
        xd    = np.linspace( xAxiso[0], xAxiso[-1], config["pic_LId"] )
        yd    = np.linspace( yAxiso[0], yAxiso[-1], config["pic_LJd"] )
        Data  = r2r.rct2rct( xp=xAxiso, yp=yAxiso, data=PEnum, xg=xd, yg=yd )
    
        
    # ------------------------------ #
    # --- [2] データ読込         --- #
    # ------------------------------ #
    Dkeys     = ["psi"]
    PEnum     = fpc.FetchPIC( job =job  , jobDir=jobDir, kstep =kstep, \
                              keys=Dkeys, config=config )
    xAxis     = PEnum["xAxis"]
    yAxis     = PEnum["yAxis"]
    flux      = PEnum["psi"  ]

    # ------------------------------ #
    # --- [3] プロット           --- #
    # ------------------------------ #
    cfs.configSettings( config=config, configType="cMap2D_def" )
    config["cmp_Nlevels"]    = 20
    config["cnt_Thick"]      = 2.0
    config["cnt_Nlevels"]    = 30
    config["cmp_ColorTable"] = "jet"
    # config["MinimalOut"]     = True
    # config["axes_x_off"]     = True
    # config["axes_y_off"]     = True
    # config["clb_sw"]         = False
    if ( Flag__PELine ):
        x1a  = yAxiso
        y1a  = np.zeros_like( x1a )
        cMap = flux
    else:
        cMap = PEnum
    cmap = clm.cMap2D( xAxis=yAxis, yAxis=xAxis, cMap=cMap, Cntr=flux,
                       FigName=OutFile, config=config, NoInstant=True )
    if ( Flag__PELine ):
        for i in range( ijDomain.shape[0] ):
            y1a[:] = xAxiso[ijDomain[i,5]]
            cmap.addPlot( xAxis=x1a, yAxis=y1a, color='Red' )
    cmap.writeFigure()
    
# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    job    = args["job"]
    jobDir = args["jobDir"]
    kstep  = args["kstep"]
    OutDir = "./png/"
    print( " [grd2cmp] ( job , kstep ) = ( {0} , {1} )".format( job, kstep ), end=""  )
    p2c    = grd2cmp( job  =job  , jobDir=jobDir, OutDir=OutDir, \
                      kstep=kstep )
    print( "\t [Done]" )
