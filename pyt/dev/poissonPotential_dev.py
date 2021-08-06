import numpy                     as np
import myStyle.cMap2D            as clm
import myStyle.sfc3D             as sfc
import mPICsee.FetchPIC          as fpc
import myStyle.LoadConfig        as lcf
import myStyle.configSettings    as cfs
import fLIB.fLIB__div2D          as div
import fLIB.fLIB__pyBiCGSTAB     as bcg

def poissonPotential( job   =None, jobDir=None, kstep=None, OutFile=None, OutDir=None, \
                      config=None, Flag__divEfigout=False, Flag__ContourOut=True, Flag__SurfaceOut=True ):
    # ------------------------- #
    # --- [1] 引数チェック  --- #
    # ------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"    .format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}png/"    .format( jobDir )
    if ( OutFile is None ): OutFile = "poisson_{0:08}.png".format( kstep  )
    OutFile = OutDir + OutFile
    
    # ------------------------- #
    # --- [2] データ読込    --- #
    # ------------------------- #
    #  -- [2-1] Data Fetch  --  #
    keys        = ["Ex","Ey","Ez","psi"]
    Data        = fpc.FetchPIC( job =job , jobDir=jobDir, kstep =kstep, \
                                keys=keys, config=config )
    xAxis       = Data["xAxis"]
    yAxis       = Data["yAxis"]
    flux        = np.ascontiguousarray( Data["psi"  ] )
    #  -- [2-2] div E       --  #
    divE        = div.div2d( u1 =Data["Ez"]   , u2 =Data["Ex"]   , coordinate="RTZ", \
                             x1g=Data["xAxis"], x2g=Data["yAxis"], difftype  ="Central__dx2" )
    #  -- [2-3] divE 描画   --  #
    if ( Flag__divEfigout ):
        FigName = OutFile.replace( "poisson", "divE" )
        clm.cMap2D( xAxis  =yAxis,   yAxis =xAxis , cMap=divE, Cntr=flux, \
                    FigName=FigName, config=config )
        
    # ------------------------- #
    # --- [3] poisson方程式 --- #
    # ------------------------- #
    #  -- [3-1] 境界条件    --  #
    divE[ 0, :] = 0.0
    divE[-1, :] = 0.0
    divE[ :, 0] = 0.0
    divE[ :,-1] = 0.0
    #  -- [3-2] poisson Inv.--  #
    rhs         = -divE
    phist       = bcg.pyBiCGSTAB( source=rhs, x1g=Data["xAxis"], x2g=Data["yAxis"], coordinate="RTZ", x2Min=0.0 )

    # ------------------------- #
    # --- [4]   プロット    --- #
    # ------------------------- #
    #  -- [4-1] Contour     --  #
    if ( Flag__ContourOut ):
        config["cmp_AutoLevel"] = False
        config["cmp_MaxMin"]    = [-0.15,0.15]
        FigName    = OutFile.replace( "poisson_", "poisson_cmp_" )
        clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis , cMap=phist, Cntr=flux, \
                    FigName=OutFile, config=config )
        
    #  -- [4-2] Surface     --  #
    if ( Flag__SurfaceOut ):
        cfs.configSettings( config=config, configType="sfc3D_def" )
        vMin, vMax = -0.15, 0.15
        FigName    = OutFile.replace( "poisson_", "poisson_sfc_" )
        sfc.sfc3D( xAxis  =yAxis  , yAxis =xAxis  , sfc =phist, Cntr=flux, \
                   elev   = +30   , azim  =+55    , vMin=vMin , vMax=vMax, \
                   FigName=FigName, config=config )
        

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    job    = args["job"]
    jobDir = args["jobDir"]
    OutDir = "./png/"
    for kstep in args["ksteps"]:
        poissonPotential( job=job, jobDir=jobDir, kstep=kstep, OutDir=OutDir )

    # # --- [3] コンフィグ --- #
    # config = lcf.LoadConfig()
    # config["AutoRange"]       = False
    # config["AutoTicks"]       = False
    # config["FigSize"]         = (5,3)
    # config["cmp_xRange"]      = [5.0,20.0]
    # config["cmp_yRange"]      = [-5.0,+5.0]
    # config["axes_x_ticks"]    = [5,10,15,20]
    # config["axes_y_ticks"]    = [-5,0,+5]
    # config["cmp_Position"]    = [0.16,0.18,+0.80,0.89]
    # config["cmp_AutoLevel"]   = True
    # config["cmp_ColorTable"]  = "jet"
    # config["clb_FontSize"]    = 12
    # config["clb_Orientation"] = "vertical"
    # config["clb_Position"]    = [0.86,0.20,0.90,0.90]
    # config["clb_Title"]       = None
    # config["MinimalOut"]      = True
    # config["cmp_MaxMin"]      = [-0.003,0.003]

    
