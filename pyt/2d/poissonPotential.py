import numpy                   as np
import myStyle.cMap2D          as clm
import myStyle.gsfc3D          as sfc
import mPICsee.FetchPIC        as fpc
import myStyle.LoadConfig      as lcf
import myStyle.configSettings  as cfs
import myConvert.where2D       as w2D
import fLIB.fLIB__div2D        as div
import fLIB.fLIB__pyBiCGSTAB   as bcg
import Filter.nxLinearFilter2D as flt
import myConvert.argmin2D      as amn

def poissonPotential( job   =None, jobDir=None, kstep=None, OutFile=None, OutDir=None, coordinate="RTZ", \
                      config=None, Flag__divEfigout=False, Flag__ContourOut=True, Flag__SurfaceOut=True ):
    # ------------------------- #
    # --- [1] 引数チェック  --- #
    # ------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}png/".format( jobDir )
    if ( OutFile is None ): OutFile = "{0}{1}poisson_{2:08}.png".format( OutDir, job, kstep )
    import mPICsee.TimeUnitConvert as tuc
    print( " kstep = {0} , time (wpi) = {1}".format( kstep, tuc.TimeUnitConvert( kstep=kstep, unit="wpi" ) ) )
    # ------------------------- #
    # --- [2] データ読込    --- #
    # ------------------------- #
    #  -- [2-1] Data Fetch  --  #
    keys        = ["Ex","Ey","Ez","psi"]
    Data        = fpc.FetchPIC( job =job , jobDir=jobDir, kstep  =kstep, \
                                keys=keys, config=config, RawData=False  )
    xAxis       = Data["xAxis"]
    yAxis       = Data["yAxis"]
    flux        = np.ascontiguousarray( Data["psi"] )
    Data        = flt.nxLinearFilter2D( Data=Data, x1Axis=xAxis, x2Axis=yAxis, keys=["Ez","Ex"],
                                        alpha=0.5, nFilter=10, coordinate="rtz" )
    #  -- [2-2] div E       --  #
    divE        = div.div2d( u1 =Data["Ez"]   , u2 =Data["Ex"]   , coordinate=coordinate, \
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
    phist       = bcg.pyBiCGSTAB( source=rhs, x1g=Data["xAxis"], x2g=Data["yAxis"], coordinate=coordinate, x2Min=0.0 )
    
    # ------------------------- #
    # --- [4]   プロット    --- #
    # ------------------------- #
    cfs.configSettings( config=config, configType="cMap2D_def" )
    vMin, vMax = -0.8, 0.8
    #  -- [4-1] Contour     --  #
    if ( Flag__ContourOut ):
        config["cmp_AutoLevel"] = False
        config["cmp_MaxMin"]    = [vMin,vMax]
        FigName    = OutFile.replace( "poisson_", "poisson_cmp_" )
        clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis , cMap=phist, Cntr=flux, \
                    FigName=OutFile, config=config )

    #  -- [4-2] Surface     --  #
    cfs.configSettings( config=config, configType="sfc3D_def" )
    config["sfc_AutoRange"] = False
    config["sfc_AutoLevel"] = True
    config["MinimalOut"]    = True
    xRange    = [-12.0,+12.0]
    yRange    = [  10.0,35.0]
    reduced   = w2D.where2D( Data=Data["psi"], x1Axis=yAxis, x2Axis=xAxis, x1Range=yRange, x2Range=xRange, Masking=False )
    flux      = reduced["Data"]
    reduced   = w2D.where2D( Data=phist, x1Axis=yAxis, x2Axis=xAxis, x1Range=yRange, x2Range=xRange, Masking=False )
    xAxis     = reduced["x2Axis"]
    yAxis     = reduced["x1Axis"]
    phist     = reduced["Data"]
    offset    = 0.2
    if ( Flag__SurfaceOut ):
        FigName    = OutFile.replace( "poisson_", "poisson_sfc_" )
        sfc.sfc3D( xAxis  =yAxis  , yAxis =xAxis , sfc=phist, Cntr=flux, cMap=phist, \
                   vMin   =vMin   , vMax  =vMax  , offset=offset, \
                   FigName=FigName, config=config )

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    OutDir = "./png/poissonPotential/"
    for kstep in args["ksteps"]:
        poissonPotential( job=args["job"], jobDir=args["jobDir"], kstep=kstep, OutDir=OutDir )
