import numpy                   as np
import mPICsee.FetchPIC        as fpc
import myStyle.LoadConfig      as lcf
import fLIB.fLIB__div2D        as div
import fLIB.fLIB__pyBiCGSTAB   as bcg
import fLIB.fLIB__grad2d       as grd
import Filter.nxLinearFilter2D as flt
    
# ======================================== #
# ===   ポアソン方程式を電場から算出   === #
# ======================================== #
def poissonPotential( job=None, jobDir=None, kstep=None, coordinate="RTZ", alpha=None, nFilter=None,
                      config=None, Flag__divEfigout=False, RowData=False ):
    # ------------------------- #
    # --- [1] 引数チェック  --- #
    # ------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( alpha   is None ): alpha   = config["cmp_LinearFilt"]
    if ( nFilter is None ): nFilter = config["cmp_nFilter"]
    # ------------------------- #
    # --- [2] データ読込    --- #
    # ------------------------- #
    #  -- [2-1] Data Fetch  --  #
    keys  = ["Ex","Ey","Ez","psi"]
    Data  = fpc.FetchPIC( job =job , jobDir=jobDir, kstep  =kstep  , \
                          keys=keys, config=config, RowData=RowData  )
    xAxis = Data["xAxis"]
    yAxis = Data["yAxis"]
    Data  = flt.nxLinearFilter2D( Data =Data , x1Axis =xAxis  , x2Axis    =yAxis, keys=["Ez","Ex"], \
                                  alpha=alpha, nFilter=nFilter, coordinate="rtz" )
    #  -- [2-2] div E       --  #
    divE  = div.div2d( u1 =Data["Ez"]   , u2 =Data["Ex"]   , coordinate=coordinate, \
                       x1g=Data["xAxis"], x2g=Data["yAxis"], difftype  ="Central__dx2" )
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
    grdp        = grd.grad2d( xg=xAxis, yg=yAxis, Data=phist )
    ret         = { "xAxis":xAxis, "yAxis":yAxis, "phist":phist, "Flux":Data["psi"] }
    ret["Estx"] = - grdp["dfdy"]
    ret["Esty"] = np.zeros_like( ret["Estx"] )
    ret["Estz"] = - grdp["dfdx"]
    ret["Eidx"] = Data["Ex"] - ret["Estx"]
    ret["Eidy"] = Data["Ey"] - ret["Esty"]
    ret["Eidz"] = Data["Ez"] - ret["Estz"]
    return( ret  )


# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    pst     = poissonPotential( job=args["job"], jobDir=args["jobDir"], kstep=args["ksteps"][0] )
    import myStyle.cMap2D as clm
    pngFile = "./png/poissonPotential_phist.png"
    clm.cMap2D( cMap=pst["phist"], FigName=pngFile )
    pngFile = "./png/poissonPotential_Estx.png"
    clm.cMap2D( cMap=pst["Estx"] , FigName=pngFile )
    pngFile = "./png/poissonPotential_Estz.png"
    clm.cMap2D( cMap=pst["Estz"] , FigName=pngFile )
    
