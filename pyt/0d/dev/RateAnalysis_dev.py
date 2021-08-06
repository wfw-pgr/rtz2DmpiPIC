import numpy                   as np
import myStyle.LoadConfig      as lcf
import myUtils.LoadConst       as lcn
import myUtils.progressBar     as pgb
import fLIB.fLIB__mpi2one      as m2o
import fLIB.fLIB__solveXOpt    as sxo
import mPICsee.calcFlux        as cfl
import myConvert.arr2dct       as a2d
import myStyle.plot1D          as pl1
import myStyle.configSettings  as cfs
import mPICsee.TimeUnitConvert as tuc
import mPICsee.FetchAxis       as fax

# =========================================================== #
# ===  レート 解析用のルーチン ( ハンドラ )               === #
# =========================================================== #
def RateAnalysis( job    =None, jobDir=None,  ksteps=None , retSize=18  , \
                  OutFile=None, OutDir=None,  config=None , pngDir =None, CnsFile=None, \
                  NewFile=True, Flag__CallCalculater=False, Flag__CallPlotter=True ):
    # ------------------------------- #
    # --- [1]   引数チェック      --- #
    # ------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( ksteps  is None ): ksteps  = np.arange( config["arg_Iter"] )*config["arg_Step"]+config["arg_Init"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}dat/".format( jobDir )
    if ( pngDir  is None ): pngDir  = "{0}png/".format( jobDir )
    if ( OutFile is None ): OutFile = OutDir + "RateAnalysis.dat"
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    # ------------------------------- #
    # --- [2] リコネクション率    --- #
    # ------------------------------- #
    if ( Flag__CallCalculater ):
        Nsteps = len( ksteps )
        if ( NewFile ): ( open( OutFile, "w" ) ).close()    # -- 新規ファイル -- #
        pgb.progressBar( iteration=0, total=Nsteps )        # -- 進捗表示     -- #
        for i,kstep in enumerate( ksteps ):
            ret = RateAnalysis_calculater( job    =job    , jobDir=jobDir, kstep=kstep, \
                                           CnsFile=CnsFile, config=config )
            pgb.progressBar( iteration=i+1, total=Nsteps )
            with open( OutFile, "ab" ) as f:
                np.savetxt( f, np.array( ret ).reshape( (1,retSize) ) )
    # ------------------------------- #
    # --- [3] 返却 /  描画        --- #
    # ------------------------------- #
    #  -- [3-1] .png 書き出し     --  #
    if ( Flag__CallPlotter ):
        RateAnalysis_plotter( datFile=OutFile, config=config, pngDir=pngDir, CnsFile=CnsFile )


# =========================================================== #
# ===  レート 解析用のルーチン ( 計算用関数 )             === #
# =========================================================== #
def RateAnalysis_calculater( job=None, jobDir=None, kstep=None, CnsFile=None, config=None, width1=5, width2=9 ):
    # -------------------------- #
    # --- [1] データ呼び出し --- #
    # -------------------------- #
    #  -- [1-1] time / rAxis --  #
    const   = lcn.LoadConst( InpFile=CnsFile )
    fstep   = float( kstep )
    axes    = fax.FetchAxis( job=job, jobDir=jobDir, kstep=kstep, CnsFile=CnsFile, unit="wce" )
    time    = axes["time"]
    zAxis   = axes["x1Axis"]
    rAxis   = axes["x2Axis"]
    #  -- [1-2] データ取得   --  #
    Data    = m2o.mpi2one( job=job, jobDir=jobDir, kstep=kstep, config=config )
    keyDict = {"Bx":0, "By":1, "Bz":2, "Bxo":6, "Byo":7, "Bzo":8, "Ex":3, "Ey":4, "ne":18, "ni":19}
    Bz      = np.ascontiguousarray( Data[ keyDict["Bz" ],:,:] )
    psiMax  = [1.0]
    flux    = cfl.calcFlux ( rg=rAxis, Bz=Bz , indexing='ji', Flag__normalize=True, fluxMax=psiMax )
    #  -- [1-3] X,O-point    --  #
    ptset   = sxo.solveXOpt( psi=flux )
    ptOXO   = ptset["ptOXO"]
    iO1, iX, iO2 = ptOXO[0,0], ptOXO[0,1], ptOXO[0,2]
    jO1, jX, jO2 = ptOXO[1,0], ptOXO[1,1], ptOXO[1,2]
    x1Xpt, x2Xpt = zAxis[iX ], rAxis[jX ]
    psiNX, psiFX = flux[jX ,iX ], flux[jX ,iX ]*psiMax[0]
    psiO1, psiO2 = flux[jO1,iO1], flux[jO2,iO2]
    psiFO        = np.max( [psiO1,psiO2] )*psiMax[0]
    # -------------------------- #
    # --- [2] Reconnection率 --- #
    # -------------------------- #
    if ( ptset["Flag__sOpt"] ):
        #  -- [2-1] Single O-point Case -- #
        Bxin  = 0.0; vAin  = 0.0; mergR  = 1.0;
        ExMR  = 0.0; EyMR  = 0.0; ErMR   = 0.0;
        ExXpt = 0.0; EyXpt = 0.0; ErXpt  = 0.0;
    else:
        #  -- [2-2]  Dual  O-point Case -- #
        #   -   dpsi/dt   -   #
        wRngj = np.array( [jX-((width1-1)/2), jX+((width1-1)/2)+1], np.int )
        wRngi = np.array( [iX-((width2-1)/2), iX+((width2-1)/2)+1], np.int )
        ExXpt = np.mean( Data[keyDict["Ex"], (wRngj[0]):(wRngj[1]), (wRngi[0]):(wRngi[1])] )
        EyXpt = np.mean( Data[keyDict["Ey"], (wRngj[0]):(wRngj[1]), (wRngi[0]):(wRngi[1])] )
        ErXpt  = np.sqrt( ExXpt**2 + EyXpt**2 )
        mergR = psiFX / psiFO        #   - Bxin / vAin -   #
        if ( mergR < 0.95 ):
            jO    = int( 0.5*( jO1 + jO2 ) )
            bx_   = Data[ keyDict["Bx"], jO,iO1:iO2 ]
            by_   = Data[ keyDict["By"], jO,iO1:iO2 ]
            bz_   = Data[ keyDict["Bz"], jO,iO1:iO2 ]
            absB  = np.sqrt( bx_**2 + by_**2 + bz_**2 )
            im    = np.argmax( absB )
            Bxin  = absB[ im ]
            wRngj = np.array( [ jO     -((width1-1)/2),  jO     +((width1-1)/2)+1], np.int )
            wRngi = np.array( [(iO1+im)-((width2-1)/2), (iO1+im)+((width2-1)/2)+1], np.int )
        else:
            iO    = int( 0.5*( iO1 + iO2 ) )
            jO    = int( 0.5*( jO1 + jO2 ) )
            wRngj = np.array( [const["LJ"]*0.1, const["LJ"]*0.9], np.int )
            bx_   = Data[ keyDict["Bx"], (wRngj[0]):(wRngj[1]), iO ]
            by_   = Data[ keyDict["By"], (wRngj[0]):(wRngj[1]), iO ]
            bz_   = Data[ keyDict["Bz"], (wRngj[0]):(wRngj[1]), iO ]
            absB  = np.sqrt( bx_**2 + by_**2 + bz_**2 )
            im    = np.argmax( absB )
            Bxin  = absB[ im ]
            wRngj = np.array( [ (wRngj[0]+im) - ((width1-1)/2), (wRngj[0]+im) + ((width1-1)/2)+1], np.int )
            wRngi = np.array( [ (iX         ) - ((width2-1)/2), (iX         ) + ((width2-1)/2)+1], np.int )
        niin  = np.mean( Data[keyDict["ni"], (wRngj[0]):(wRngj[1]), (wRngi[0]):(wRngi[1])] )
        vAin  = Bxin / ( const["wpewce"]*np.sqrt( const["mr"]*niin ) )
        #   - Reconnection Rate -   #
        import myBasicAlgs.robustInv as inv
        ExMR  = ExXpt * inv.robustInv( vAin * Bxin )
        EyMR  = EyXpt * inv.robustInv( vAin * Bxin )
        ErMR  = ErXpt * inv.robustInv( vAin * Bxin )
    return( [ fstep, time , psiFX, psiFO, psiNX, psiO1, psiO2, x1Xpt, x2Xpt, vAin, Bxin, \
              ExXpt, EyXpt, ErXpt, mergR, ExMR , EyMR , ErMR  ] )


# --- .png 書き出し --- #
def RateAnalysis_plotter( datFile=None, config=None, pngDir=None, CnsFile=None ):
    # ---------------------------------------- #
    # --- [1] プロット準備                  -- #
    # ---------------------------------------- #
    #  -- [1-1] データ読込                 --  #
    Data  = np.transpose( np.loadtxt( datFile ) )
    keys  = [ "fstep", "time" , "psiFX", "psiFO", "psiNX", "psiO1", "psiO2", \
              "x1Xpt", "x2Xpt", "vAin" , "Bxin" , "ExXpt", "EyXpt", "ErXpt", \
              "mergR", "ExMR" , "EyMR" , "ErMR" ]
    Data  = a2d.arr2dct( Data=Data, keys=keys )
    ptime = tuc.TimeUnitConvert( Data=Data["time"], unit="wpi", CnsFile=CnsFile )
    MRate = ( Data["mergR"] - Data["mergR"][0] ) / ( Data["mergR"][-1] - Data["mergR"][0] )
    yRate = - Data["EyMR" ]
    #  -- [1-1] コンフィグ                 --  #
    cfs.configSettings( config=config, configType="plot1D_lateral" )
    cfs.configSettings( config=config, configType="plot1D_stack" )
    config["plt_xAutoRange"]  = False
    config["plt_xRange"]      = [0.,120.]
    config["plt_yAutoRange"]  = False

    # ---------------------------------------- #
    # --- [2] プロット                     --- #
    # ---------------------------------------- #
    #  --      Merging Rate -- #
    cfs.configSettings( config=config, configType="FilterClear" )
    OutFile = pngDir + "t-MergingRate.png"
    config["yMajor_Nticks"]   = 7
    config["plt_yRange"]      = [0.0,1.2]
    config["plt_GaussFilt"]   = 1.2
    config["plt_color"]       = "green"
    pl1.plot1D( xAxis=ptime, yAxis=MRate, config=config, FigName=OutFile )
    #  -- Reconnection Rate -- #
    cfs.configSettings( config=config, configType="FilterClear" )
    OutFile = pngDir + "t-yMRRate.png"
    config["yMajor_Nticks"]   = 5
    config["plt_yRange"]      = [-0.1,0.3]
    config["plt_color"]       = "crimson"
    config["cursor_y"]        = [0.0]
    config["cursor_style"]    = "-"
    config["cursor_width"]    = 0.8
    config["cursor_color"]    = "dimgrey"
    pl1.plot1D( xAxis=ptime, yAxis=yRate, config=config, FigName=OutFile )
    

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    pngDir  = "./png/"
    print( "[RateAnalysis] {0} is under processed...".format( args["job"] ), end="\n" )
    Output  = RateAnalysis( job=args["job"], pngDir=pngDir, ksteps=args["ksteps"] )
    print( "\t [Done]" )
