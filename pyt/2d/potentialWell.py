import numpy                   as np
import mPICsee.FetchPIC        as fpc
import myStyle.LoadConfig      as lcf
import myConvert.where2D       as w2D
import fLIB.fLIB__div2D        as div
import fLIB.fLIB__pyBiCGSTAB   as bcg
import Filter.nxLinearFilter2D as flt
import fLIB.fLIB__solveXOpt    as sxo
import mPICsee.TimeUnitConvert as tuc

def potentialWell( job   =None, jobDir=None, kstep=None, OutFile=None, OutDir=None, coordinate="RTZ", \
                   config=None, Flag__1DplotOut=True   , avgRange=[-1.0,1.0], width1=2, width2=5, mRate0=None ):
    
    # ------------------------- #
    # --- [1] 引数チェック  --- #
    # ------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "out/"
    if ( OutFile is None ): OutFile = "{0}{1}pWell{2:08}.dat".format( OutDir, job, kstep )
    
    # ------------------------- #
    # --- [2] データ読込    --- #
    # ------------------------- #
    #  -- [2-1] Data Fetch  --  #
    keys  = ["Ex","Ey","Ez","psi"]
    Data  = fpc.FetchPIC( job =job , jobDir=jobDir, kstep  =kstep, \
                          keys=keys, config=config, RawData=False   )
    xAxis = Data["xAxis"]
    yAxis = Data["yAxis"]
    Data  = flt.nxLinearFilter2D( Data=Data, x1Axis=xAxis, x2Axis=yAxis, keys=["Ez","Ex"],
                                  alpha=0.2, nFilter=5, coordinate="rtz" )
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
    
    # ---------------------------------------- #
    # --- [4] 抜き出し                     --- #
    # ---------------------------------------- #
    # -- time    -- #
    fstep,time  = float( kstep ), tuc.TimeUnitConvert( job=job, jobDir=jobDir, kstep=kstep )
    # -- 1D plot -- #
    reduced     = w2D.where2D( Data=phist, x2Axis=xAxis, x2Range=avgRange )
    phi1D       = np.ravel( np.mean( reduced["Data"], axis=1 )  )
    if ( Flag__1DplotOut ):
        wData       = np.zeros( (yAxis.size,2) )
        wData[:,0]  = yAxis
        wData[:,1]  = phi1D
        np.savetxt( OutFile, wData )
    # --   V@X   -- #
    ptOXO       = ( sxo.solveXOpt( psi=Data["psi"] ) )["ptOXO"]
    iX,  jX     = ptOXO[0,1], ptOXO[1,1]
    iX1, iX2    = iX-width1, iX+width1
    jX1, jX2    = jX-width2, jX+width2
    meanV       = np.mean( phist[jX1:(jX2+1),iX1:(iX2+1)] )
    rX,  zX     = yAxis[jX], xAxis[iX]
    mRate       = ( Data["psi"] )[jX,iX] / np.max( [ ( Data["psi"] )[ptOXO[1,0],ptOXO[0,0]], ( Data["psi"] )[ptOXO[1,2],ptOXO[0,2]] ] )
    if ( mRate0[0] < 0.0 ): mRate0[0] = mRate
    mRate       = mRate - mRate0[0]
    # --   minV  -- #
    jMin, iMin  = np.argmin( phi1D ), int( ( phist.shape[1] )*0.5 )
    VMin        = phi1D[jMin]
    rMin, zMin  = yAxis[jMin], xAxis[iMin]
    # --   Eeff  -- #
    Eeff  = ( meanV - VMin ) / ( rX - rMin )
    # --   出力  -- #
    wData       = np.array( [ fstep, time, mRate, iX, jX, rX, zX, meanV, iMin, jMin, rMin, zMin, VMin, meanV-VMin, Eeff ] )
    return( wData )
    
# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    OutDir  = "./potentialWell/"
    datFile = "potentialDiff.dat"
    NewFile = True
    mRate0  = [-1.0]
    if ( NewFile ):
        ( open( datFile, "w" ) ).close()
    for kstep in args["ksteps"]:
        ret = potentialWell( job=args["job"], jobDir=args["jobDir"], kstep=kstep, OutDir=OutDir, mRate0=mRate0 )
        with open( datFile, "ab" ) as f:
            np.savetxt( f, ret.reshape( 1,-1 ) )
