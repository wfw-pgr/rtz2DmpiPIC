import numpy                   as np
import myStyle.cMap2D          as clm
import mPICsee.FetchPIC        as fpc
import myStyle.LoadConfig      as lcf
import myStyle.configSettings  as cfs
import fLIB.fLIB__div2D        as div
import fLIB.fLIB__pyBiCGSTAB   as bcg
import fLIB.fLIB__grad2d       as grd
import Filter.nxLinearFilter2D as flt

# ======================================== #
# ===  Static / Inductive Field を描く === #
# ======================================== #
def StaticOrInductive( job  =None, jobDir =None, kstep =None, OutFile=None, OutDir=None, coordinate="RTZ", \
                       alpha=0.1 , nFilter=10, config=None, Flag__divEfigout=False, Flag__ContourOut=True ):
    # ------------------------- #
    # --- [1] 引数チェック  --- #
    # ------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}png/".format( jobDir )
    if ( OutFile is None ): OutFile = "{0}out_{1}_{2:08}.png".format( OutDir, job, kstep )
    import mPICsee.TimeUnitConvert as tuc
    print( " kstep = {0} , time (wpi) = {1}".format( kstep, tuc.TimeUnitConvert( kstep=kstep, unit="wpi" ) ) )
    # ------------------------- #
    # --- [2] データ読込    --- #
    # ------------------------- #
    #  -- [2-1] Data Fetch  --  #
    keys  = ["Ex","Ey","Ez","psi"]
    Data  = fpc.FetchPIC( job =job , jobDir=jobDir, kstep  =kstep, \
                          keys=keys, config=config, RawData=False  )
    xAxis = Data["xAxis"]
    yAxis = Data["yAxis"]
    flux  = np.ascontiguousarray( Data["psi"] )
    Data  = flt.nxLinearFilter2D( Data=Data  , x1Axis =xAxis  , x2Axis=yAxis, keys=["Ez","Ex"],
                                  alpha=alpha, nFilter=nFilter, coordinate="rtz" )
    #  -- [2-2] div E       --  #
    divE  = div.div2d( u1 =Data["Ez"]   , u2 =Data["Ex"]   , coordinate=coordinate, \
                       x1g=Data["xAxis"], x2g=Data["yAxis"], difftype  ="Central__dx2" )
    #  -- [2-3] divE 描画   --  #
    if ( Flag__divEfigout ):
        FigName = OutFile.replace( "out", "divE" )
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
    rhs         = - divE
    phist       = bcg.pyBiCGSTAB( source=rhs, x1g=Data["xAxis"], x2g=Data["yAxis"], \
                                  coordinate=coordinate, x2Min=0.0 )
    #  -- [3-3] 電場再計算  --  #
    grdp        = grd.grad2d( xg=xAxis, yg=yAxis, Data=phist )
    Est, Eid    = {}, {}
    Est["Ex"]   = - grdp["dfdy"]
    Est["Ez"]   = - grdp["dfdx"]
    Eid["Ex"]   = Data["Ex"] - Est["Ex"]
    Eid["Ez"]   = Data["Ez"] - Est["Ez"]
    
    # ------------------------- #
    # --- [4]   プロット    --- #
    # ------------------------- #
    cfs.configSettings( config=config, configType="cMap2D_thesis" )
    cfs.configSettings( config=config, configType="vector2D_def"  )
    cfs.configSettings( config=config, configType="NoAxis"        )
    config["vec_scale"]      = 0.03
    config["vec_color"]      = "aqua"
    config["cmp_AutoLevel"]  = True
    config["cmp_ColorTable"] = "plasma"
    config["clb_sw"]         = False
    # -- ポテンシャル -- #
    FigName    = OutFile.replace( "out_", "phis_" )
    clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis , cMap=phist, Cntr=flux, \
                xvec   =Data["Ex"], yvec   =Data["Ez"], \
                FigName=FigName, config=config )
    # --  静電場  -- #
    config["cmp_AutoLevel"] = False
    config["cmp_MaxMin"]    = [-0.2,0.2]
    FigName    = OutFile.replace( "out_", "Estz_" )
    clm.cMap2D( xAxis  =yAxis    , yAxis =xAxis    , cMap=Est["Ez"], Cntr=flux, \
                xvec   =Est["Ex"], yvec  =Est["Ez"], \
                FigName=FigName  , config=config     )
    FigName    = OutFile.replace( "out_", "Estx_" )
    clm.cMap2D( xAxis  =yAxis    , yAxis =xAxis    , cMap=Est["Ex"], Cntr=flux, \
                xvec   =Est["Ex"], yvec  =Est["Ez"], \
                FigName=FigName  , config=config     )
    # -- 誘導電場 -- #
    config["cmp_MaxMin"]    = [-0.1,0.1]
    FigName    = OutFile.replace( "out_", "Eidz_" )
    clm.cMap2D( xAxis  =yAxis    , yAxis =xAxis    , cMap=Eid["Ez"], Cntr=flux, \
                xvec   =Eid["Ex"], yvec  =Eid["Ez"], \
                FigName=FigName  , config=config     )
    FigName    = OutFile.replace( "out_", "Eidx_" )
    clm.cMap2D( xAxis  =yAxis    , yAxis =xAxis    , cMap=Eid["Ex"], Cntr=flux, \
                xvec   =Eid["Ex"], yvec  =Eid["Ez"], \
                FigName=FigName  , config=config     )
    # -- 合計電場 -- #
    config["cmp_MaxMin"]    = [-0.2,0.2]
    FigName    = OutFile.replace( "out_", "Ettz_" )
    clm.cMap2D( xAxis  =yAxis    , yAxis =xAxis    , cMap=Data["Ez"], Cntr=flux, \
                xvec   =Eid["Ex"], yvec  =Eid["Ez"], \
                FigName=FigName  , config=config     )
    FigName    = OutFile.replace( "out_", "Ettx_" )
    clm.cMap2D( xAxis  =yAxis    , yAxis =xAxis    , cMap=Data["Ex"], Cntr=flux, \
                xvec   =Eid["Ex"], yvec  =Eid["Ez"], \
                FigName=FigName  , config=config     )

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    OutDir = "./png/StaticOrInductive/"
    for kstep in args["ksteps"]:
        StaticOrInductive( job=args["job"], jobDir=args["jobDir"], kstep=kstep, OutDir=OutDir )
