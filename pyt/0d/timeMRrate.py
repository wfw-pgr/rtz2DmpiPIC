import numpy                as np
import myStyle.LoadConfig   as lcf
import myUtils.LoadConst    as lcn
import myUtils.progressBar  as pgb
import fLIB.fLIB__mpi2one   as m2o
import fLIB.fLIB__solveXOpt as sxo
import mPICsee.calcFlux     as cfl
import myConvert.arr2dct    as a2d

# =========================================================== #
# ===  レート 解析用のルーチン ( ハンドラ )               === #
# =========================================================== #
def timeMRrate( job    =None, jobDir=None, ksteps=None, \
                OutFile=None, OutDir=None, config=None, pngDir=None, CnsFile=None, \
                NewFile=True, Flag__CallCalculater=True, Flag__CallPlotter=True ):
    # ------------------------------- #
    # --- [1]   引数チェック      --- #
    # ------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( ksteps  is None ): ksteps  = np.arange( config["arg_Iter"] )*config["arg_Step"]+config["arg_Init"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}dat/".format( jobDir )
    if ( pngDir  is None ): pngDir  = "{0}png/".format( jobDir )
    if ( OutFile is None ): OutFile = OutDir + "timeMRrate.dat"
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    # ------------------------------- #
    # --- [2] リコネクション率    --- #
    # ------------------------------- #
    if ( Flag__CallCalculater ):
        Nsteps = len( ksteps )
        if ( NewFile ): ( open( OutFile, "w" ) ).close()    # -- 新規ファイル -- #
        pgb.progressBar( iteration=0, total=Nsteps )        # -- 進捗表示     -- #
        for i,kstep in enumerate( ksteps ):
            ret = timeMRrate_calculater( job    =job    , jobDir=jobDir, kstep=kstep, \
                                         CnsFile=CnsFile, config=config )
            pgb.progressBar( iteration=i+1, total=Nsteps )
            with open( OutFile, "ab" ) as f:
                np.savetxt( f, np.array( ret ).reshape( (1,13) ) )
    
    # ------------------------------- #
    # --- [3] 返却 /  描画        --- #
    # ------------------------------- #
    #  -- [3-1] .png 書き出し     --  #
    if ( Flag__CallPlotter ):
        timeMRrate_plotter( datFile=OutFile, config=config, OutDir=pngDir )


def timeMRrate_calculater( job=None, jobDir=None, kstep=None, CnsFile=None, config=None, width1=5, width2=9 ):
    # -------------------------- #
    # --- [1] データ呼び出し --- #
    # -------------------------- #
    #  -- [1-1] time / rAxis --  #
    const   = lcn.LoadConst( InpFile=CnsFile )
    dtInv   =            1.0 / const["dt"]
    time    = float( kstep ) * const["dt"]
    rAxis   = np.linspace( const["x2Min"], const["x2Min"]+const["x2Leng"], const["LJ"]  )
    #  -- [1-2] データ取得   --  #
    Data    = m2o.mpi2one( job=job, jobDir=jobDir, kstep=kstep, config=config )
    keyDict = {"Bx":0, "By":1, "Bz":2, "Bxo":6, "Byo":7, "Bzo":8, "Ex":3, "Ey":4, "ne":18, "ni":19}
    Bz      = np.ascontiguousarray( Data[ keyDict["Bz" ],:,:] )
    Bzo     = np.ascontiguousarray( Data[ keyDict["Bzo"],:,:] )
    psiMax1 = [1.0]
    psiMax2 = [1.0]
    flux1   = cfl.calcFlux ( rg=rAxis, Bz=Bz , indexing='ji', Flag__normalize=True, fluxMax=psiMax1 )
    flux2   = cfl.calcFlux ( rg=rAxis, Bz=Bzo, indexing='ji', Flag__normalize=True, fluxMax=psiMax2 )
    #  -- [1-3] X,O-point    --  #
    set1    = sxo.solveXOpt( psi=flux1 )
    set2    = sxo.solveXOpt( psi=flux2 )
    ptOXO1  = set1["ptOXO"]
    ptOXO2  = set2["ptOXO"]
    iO1_1, iX_1, iO2_1  = ptOXO1[0,0], ptOXO1[0,1], ptOXO1[0,2]
    jO1_1, jX_1, jO2_1  = ptOXO1[1,0], ptOXO1[1,1], ptOXO1[1,2]
    iO1_2, iX_2, iO2_2  = ptOXO2[0,0], ptOXO2[0,1], ptOXO2[0,2]
    jO1_2, jX_2, jO2_2  = ptOXO2[1,0], ptOXO2[1,1], ptOXO2[1,2]
    # -------------------------- #
    # --- [2] Reconnection率 --- #
    # -------------------------- #
    if ( set1["Flag__sOpt"] or set2["Flag__sOpt"] ):
        #  -- [2-1] Single O-point Case -- #
        Bxin   = 0.0; vAin  = 0.0
        dpdt1  = 0.0; E_MR  = 0.0
        mrgR   = 1.0;
        dpdt2  = 0.0; EXpt  = 0.0
    else:
        #  -- [2-2]  Dual  O-point Case -- #
        #   -   dpsi/dt   -   #
        rXpt  = rAxis[ jX_1 ]
        psiX1 = flux1[ jX_1 , iX_1 ]
        psiX2 = flux2[ jX_2 , iX_2 ]
        psiO1 = np.max([ flux1[ jO1_1, iO1_1], flux1[ jO2_1, iO2_1] ])
        psiO2 = np.max([ flux2[ jO1_2, iO1_2], flux2[ jO2_2, iO2_2] ])
        psiX  = psiX1*psiMax1[0]
        psiO  = psiO1*psiMax1[0]
        dpdt1 = ( ( psiX2-psiO2 )*psiMax2[0] - ( ( psiX1-psiO1 )*psiMax1[0] ) ) * dtInv
        dpdt2 = ( psiX1*psiMax1[0] - psiX2*psiMax2[0] )*dtInv
        wRngj = np.array( [jX_1-((width1-1)/2), jX_1+((width1-1)/2)+1], np.int )
        wRngi = np.array( [iX_1-((width2-1)/2), iX_1+((width2-1)/2)+1], np.int )
        ExXpt = np.mean( ( Data[keyDict["Ex"]] )[(wRngj[0]):(wRngj[1]),(wRngi[0]):(wRngi[1])] )
        EyXpt = np.mean( ( Data[keyDict["Ey"]] )[(wRngj[0]):(wRngj[1]),(wRngi[0]):(wRngi[1])] )
        mrgR  =  ( psiX1+psiX2 ) / ( psiO1+psiO2 )
        #   - Bxin / vAin -   #
        jO_1  = int( 0.5*( jO1_1 + jO2_1 ) )
        bx_   =  ( 0.5*( Data[keyDict["Bx"]] + Data[keyDict["Bxo"]] ) )[ jO_1,iO1_1:iO2_1 ]
        by_   =  ( 0.5*( Data[keyDict["By"]] + Data[keyDict["Byo"]] ) )[ jO_1,iO1_1:iO2_1 ]
        # ni_   =  (       Data[keyDict["ni"]]                          )[ jO_1,iO1_1:iO2_1 ]
        absB  = np.sqrt( bx_**2 + by_**2 )
        im    = np.argmax( absB )
        Bxin  = absB[ im ]
        wRngj = np.array( [ jO_1     -((width1-1)/2),  jO_1     +((width1-1)/2)+1], np.int )
        wRngi = np.array( [(iO1_1+im)-((width1-1)/2), (iO1_1+im)+((width1-1)/2)+1], np.int )
        niin  = np.mean( ( Data[keyDict["ni"]] )[(wRngj[0]):(wRngj[1]),(wRngi[0]):(wRngi[1])] )
        vAin  = Bxin / ( const["wpewce"]*np.sqrt( const["mr"]*niin ) )
        #   - Reconnection Rate -   #
        E_MR  =  dpdt1 / ( ( 2.0*np.pi*rXpt )*( vAin*Bxin ) )
        EXpt  =  dpdt2 / ( ( 2.0*np.pi*rXpt )*( vAin*Bxin ) )
        return( [ time, vAin, Bxin, psiX, psiO, rXpt, dpdt1, dpdt2, ExXpt, EyXpt, mrgR, EXpt, E_MR  ] )


# --- .png 書き出し --- #
def timeMRrate_plotter( datFile=None, config=None, OutDir=None ):
    # ------------------------- #
    # --- [1] プロット準備   -- #
    # ------------------------- #
    #  -- [1-1] データ読込  --  #
    Data  = np.transpose( np.loadtxt( datFile ) )
    Data  = removeIrregular( Data=Data )
    keys  = [ "time", "vAin", "Bxin", "psiX", "psiO", "rXpt", "dpdt1", "dpdt2", "ExXpt", "EyXpt", "mrgR", "EXpt", "E_MR"]
    label = [ "Time", "$v_{Ain}$", "B_{in}", "$\Psi_{X}$", "$\Psi_{O}$", "$r_{X}$", "$d\psi /dt$", "$d\psi /dt$",
              "\eta_{Merg}", "$E_{x}$", "$E_{y}$", "$E_{Xpt}$", "$E_R$" ]
    Data  = a2d.arr2dct( Data=Data , keys=keys )
    label = a2d.arr2dct( Data=label, keys=keys )
    #  -- [1-2] コンフィグ   -- #
    import myStyle.plot1D         as pl1
    import myStyle.configSettings as cfs
    cfs.configSettings( config=config, configType="plot1D_def" )
    config["xTitle"]          = "Time"
    config["yTitle"]          = ""
    config["plt_LegFontSize"] = 10
    config["plt_LegNColumn"]  = 1
    config["plt_xAutoRange"]  = False
    config["plt_xRange"]      = [0.,200.]

    # ----------------------- #
    # --- [2] プロット    --- #
    # ----------------------- #
    #  -- [2-1] vAin, Bxin, dpdt  -- #
    OutFile = OutDir + "timeMergingRate.png"
    pl1.plot1D( xAxis=Data["time"], yAxis=Data["mrgR"], label=label["mrgR"], config=config, FigName=OutFile )
    OutFile = OutDir + "timeEXpt.png"
    pl1.plot1D( xAxis=Data["time"], yAxis=Data["EXpt"], label=label["EXpt"], config=config, FigName=OutFile )
    OutFile = OutDir + "timedpdt1.png"
    pl1.plot1D( xAxis=Data["time"], yAxis=Data["dpdt1"], label=label["dpdt1"], config=config, FigName=OutFile )
    OutFile = OutDir + "timedpdt2.png"
    pl1.plot1D( xAxis=Data["time"], yAxis=Data["dpdt2"], label=label["dpdt2"], config=config, FigName=OutFile )
    OutFile = OutDir + "timeExXpt.png"
    pl1.plot1D( xAxis=Data["time"], yAxis=Data["ExXpt"], label=label["ExXpt"], config=config, FigName=OutFile )
    OutFile = OutDir + "timeEyXpt.png"
    pl1.plot1D( xAxis=Data["time"], yAxis=Data["EyXpt"], label=label["EyXpt"], config=config, FigName=OutFile )
    #  -- [2-2] リコネクション率  -- #
    config["plt_yAutoRange"] = False
    config["plt_yRange"]     = [0.,10.0]
    OutFile = OutDir + "timeMRrate.png"
    pl1.plot1D( xAxis=Data["time"], yAxis=Data["E_MR"], label=label["E_MR"], config=config, FigName=OutFile )
    #  -- [2-3] RecAvgRate -- #
    import fLIB.fLIB__derivative as bbn
    vAin = Data["vAin"]
    Bxin = Data["Bxin"]
    vAin[np.where( vAin<0.05 ) ] = 0.05
    Bxin[np.where( Bxin<1.00 ) ] = 1.00
    coef = 1.0 / ( 2.0*np.pi*Data["rXpt"]*vAin*Bxin )
    coef_ = np.zeros( coef.size-1 )
    vBin = vAin * Bxin
    for i in range( coef_.size ): coef_[i] = ( coef[i] + coef[i+1] )*0.5
    time_, Erec = bbn.myDerivative( xAxis=Data["time"], yAxis=Data["psiX"] )
    Erec_ = coef_ * Erec
    config["plt_yAutoRange"]  = True
    OutFile = OutDir + "timeErec.png"
    pl1.plot1D( xAxis=time_, yAxis=Erec, label=label["E_MR"], config=config , FigName=OutFile )
    OutFile = OutDir + "timeVBin.png"
    pl1.plot1D( xAxis=Data["time"], yAxis=vBin, config=config , FigName=OutFile )
    OutFile = OutDir + "timeMRrateAveraged.png"
    pl1.plot1D( xAxis=time_, yAxis=Erec_, label=label["E_MR"], config=config, FigName=OutFile )

    config["plt_yAutoRange"] = False
    config["plt_yRange"]     = [-0.2,0.2]
    ExXpt = - Data["ExXpt"] / ( vAin*Bxin )
    EyXpt = - Data["EyXpt"] / ( vAin*Bxin )
    OutFile = OutDir + "timeMRrate_x.png"
    pl1.plot1D( xAxis=Data["time"], yAxis=ExXpt, label=label["ExXpt"], config=config, FigName=OutFile )
    OutFile = OutDir + "timeMRrate_y.png"
    pl1.plot1D( xAxis=Data["time"], yAxis=EyXpt, label=label["EyXpt"], config=config, FigName=OutFile )
    ExyXpt  = np.sqrt( ExXpt**2 + EyXpt**2 )
    OutFile = OutDir + "timeMRrate_xy.png"
    pl1.plot1D( xAxis=Data["time"], yAxis=ExyXpt, label=label["EyXpt"], config=config, FigName=OutFile )
    config["plt_yAutoRange"] = True
    OutFile = OutDir + "vAin.png"
    pl1.plot1D( xAxis=Data["time"], yAxis=vAin, label=label["vAin"], config=config, FigName=OutFile )
    OutFile = OutDir + "Bxin.png"
    pl1.plot1D( xAxis=Data["time"], yAxis=Bxin, label=label["Bxin"], config=config, FigName=OutFile )
    
    

def removeIrregular( Data=None ):
    if ( Data is None ): sys.exit( "[timeMRrate] None Data" )
    return( Data )

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    OutDir  = None
    print( "[timeMRrate] {0} is under processed...".format( args["job"] ), end="\n" )
    Output  = timeMRrate( job=args["job"], OutDir=OutDir, ksteps=args["ksteps"] )
    print( "\t [Done]" )
