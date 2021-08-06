import mPICsee.FetchPIC        as fpc
import myStyle.LoadConfig      as lcf
import fLIB.fLIB__div2D        as div
import fLIB.fLIB__pyBiCGSTAB   as bcg
import fLIB.fLIB__grad2d       as grd

# ======================================== #
# ===  Static / Inductive Field を描く === #
# ======================================== #
def EstidDecompose( job    =None , jobDir    =None , kstep =None, \
                    RawData=False, coordinate="RTZ", config=None ):
    # ------------------------- #
    # --- [1] 引数チェック  --- #
    # ------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    # ------------------------- #
    # --- [2] データ読込    --- #
    # ------------------------- #
    #  -- [2-1] Data Fetch  --  #
    keys  = ["Ex","Ey","Ez","psi"]
    Data  = fpc.FetchPIC( job =job , jobDir    =jobDir, kstep =kstep, \
                          keys=keys, FilterKeys=keys  , config=config, RawData=RawData  )
    xAxis = Data["xAxis"]
    yAxis = Data["yAxis"]
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
    rhs         = - divE
    phist       = bcg.pyBiCGSTAB( source=rhs, x1g=Data["xAxis"], x2g=Data["yAxis"], \
                                  coordinate=coordinate, x2Min=0.0 )
    #  -- [3-3] 電場再計算  --  #
    grdp        = grd.grad2d( xg=xAxis, yg=yAxis, Data=phist )
    ret         = {}
    ret["Estx"] = - grdp["dfdy"]
    ret["Estz"] = - grdp["dfdx"]
    ret["Eidx"] = Data["Ex"] - ret["Estx"]
    ret["Eidz"] = Data["Ez"] - ret["Estz"]
    return( ret )

# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    ret     = EstidDecompose( job=args["job"], jobDir=args["jobDir"] )
    print( ret.keys() )
