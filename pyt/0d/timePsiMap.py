import numpy                as np
import myStyle.LoadConfig   as lcf
import myUtils.LoadConst    as lcn
import myUtils.progressBar  as pgb
import fLIB.fLIB__mpi2one   as m2o
import fLIB.fLIB__solveXOpt as sxo
import mPICsee.calcFlux     as cfl
import myConvert.arr2dct    as a2d

def timePsiMap( job    =None, jobDir=None, ksteps=None , \
                OutFile=None, OutDir=None, config=None , pngDir=None, CnsFile=None, \
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
    if ( OutFile is None ): OutFile = OutDir + "timePsiMap.dat"
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    # ------------------------------- #
    # --- [2] Psi Map を0次元化   --- #
    # ------------------------------- #
    if ( Flag__CallCalculater ):
        psiX0  = []
        Nsteps = len( ksteps )
        if ( NewFile ): ( open( OutFile, "w" ) ).close()    # -- 新規ファイル -- #
        pgb.progressBar( iteration=0, total=Nsteps )        # -- 進捗表示     -- #
        for i,kstep in enumerate( ksteps ):
            ret = timePsiMap_calculater( job    =job    , jobDir=jobDir, kstep=kstep, \
                                         CnsFile=CnsFile, config=config, psiX0=psiX0 )
            pgb.progressBar( iteration=i+1, total=Nsteps )
            with open( OutFile, "ab" ) as f:
                np.savetxt( f, np.array( ret ).reshape( (1,12) ) )
    
    # ------------------------------- #
    # --- [3] 返却 /  描画        --- #
    # ------------------------------- #
    #  -- [3-1] .png 書き出し     --  #
    if ( Flag__CallPlotter ):
        timePsiMap_plotter( datFile=OutFile, config=config, OutDir=pngDir )


def timePsiMap_calculater( job=None, jobDir=None, kstep=None, CnsFile=None, config=None, psiX0=None ):
    # -------------------------- #
    # --- [1] データ呼び出し --- #
    # -------------------------- #
    #  -- [1-1] time / rAxis --  #
    const   = lcn.LoadConst( InpFile=CnsFile )
    time    = kstep * const["dt"]
    zAxis   = np.linspace( -0.5*const["x1Leng"], 0.5*const["x1Leng"], const["LI"]  )
    rAxis   = np.linspace( const["x2Min"], const["x2Min"]+const["x2Leng"], const["LJ"]  )
    #  -- [1-2] データ取得   --  #
    Data    = m2o.mpi2one( job=job, jobDir=jobDir, kstep=kstep, config=config )
    keyDict = { "Bz":2 }
    Bz      = np.ascontiguousarray( Data[ keyDict["Bz"],:,:] )
    psiMax  = 1.0
    flux    = cfl.calcFlux ( rg=rAxis, Bz=Bz , indexing='ji', Flag__normalize=True, fluxMax=psiMax )
    #  -- [1-3] X,O-point    --  #
    setPsi  = sxo.solveXOpt( psi=flux )
    ptOXO   = setPsi["ptOXO"]
    iO1, iX, iO2  = ptOXO[0,0], ptOXO[0,1], ptOXO[0,2]
    jO1, jX, jO2  = ptOXO[1,0], ptOXO[1,1], ptOXO[1,2]
    # -------------------------- #
    # --- [2] Psi のマップ   --- #
    # -------------------------- #
    #   -   psi    -   #
    psiFX = flux[ jX , iX ]
    psiO1 = flux[ jO1, iO1]
    psiO2 = flux[ jO2, iO2]
    psiNX = psiFX / ( 0.5*(psiO1 + psiO2) )
    if ( len( psiX0 ) == 0 ): psiX0.append( psiNX )
    mrgR  = ( psiNX - psiX0[0] ) / ( 1.0 - psiX0[0] )
    #   - Position -   #
    rOpt1 = rAxis[jO1]
    rXpt  = rAxis[jX ]
    rOpt2 = rAxis[jO2]
    zOpt1 = zAxis[iO1]
    zXpt  = zAxis[iX ]
    zOpt2 = zAxis[iO2]
    return( [ time, psiO1, psiO2, psiFX, psiNX, rOpt1, zOpt1, rOpt2, zOpt2, rXpt, zXpt, mrgR ] )


# --- .png 書き出し --- #
def timePsiMap_plotter( datFile=None, config=None, OutDir=None ):
    # ------------------------- #
    # --- [1] プロット準備   -- #
    # ------------------------- #
    #  -- [1-1] データ読込  --  #
    Data  = np.transpose( np.loadtxt( datFile ) )
    keys  = [ "time", "psiO1", "psiO2", "psiFX", "psiNX",
              "rOpt1", "zOpt1", "rOpt2", "zOpt2", "rXpt", "zXpt", "mrgR" ]
    label = [ "Time", "$\Psi_{O1}$", "$\Psi_{O2}$", "$\Psi_{X}$", "$\Psi_{Xnrm}$",
              "$r_{O1}$", "$z_{O1}$", "$r_{O2}$", "$z_{O2}$", "$r_{X}$", "$z_{X}$", "$\eta_{Merg}$" ]
    Data  = a2d.arr2dct( Data=Data , keys=keys )
    label = a2d.arr2dct( Data=label, keys=keys )
    #  -- [1-2] コンフィグ   -- #
    import myStyle.plot1D         as pl1
    import myStyle.configSettings as cfs
    cfs.configSettings( config=config, configType="plot1D_def" )
    config["xTitle"]          = "Time"
    config["yTitle"]          = ""
    config["plt_xAutoRange"]  = False
    config["plt_yAutoRange"]  = False
    config["plt_xRange"]      = [0.0,180.0]
    config["plt_LegFontSize"] = 16
    config["plt_LegNColumn"]  = 3

    # ----------------------- #
    # --- [2] プロット    --- #
    # ----------------------- #
    #  -- [2-1] psi Value  -- #
    OutFile = OutDir + "timePsiMap_psi.png"
    config["plt_yRange"]      = [0.0,1.2]
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["psiO1"], label=label["psiO1"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["psiO2"], label=label["psiO2"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["psiFX"], label=label["psiFX"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    #  -- [2-2] rOpt, rXpt -- #
    OutFile = OutDir + "timePsiMap_rOXpt.png"
    config["plt_yRange"]      = [2.0,6.0]
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["rOpt1"], label=label["rOpt1"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["rOpt2"], label=label["rOpt2"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["rXpt" ], label=label["rXpt" ] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    #  -- [2-3] rOpt, rXpt -- #
    OutFile = OutDir + "timePsiMap_zOXpt.png"
    config["plt_yRange"]      = [-4.0,4.0]
    fig     = pl1.plot1D( config=config, FigName=OutFile )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["zOpt1"], label=label["zOpt1"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["zOpt2"], label=label["zOpt2"] )
    fig.addPlot( xAxis=Data["time"], yAxis=Data["zXpt" ], label=label["zXpt" ] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    #  -- [2-4] mrgR -- #
    OutFile = OutDir + "timePsiMap_mrgR.png"
    config["yTitle"]          = "$\eta_{Merg}$"
    config["plt_yRange"]      = [0.0,1.2]
    fig     = pl1.plot1D( xAxis =Data["time"], yAxis  =Data["mrgR" ], \
                          config=config      , FigName=OutFile )


# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    OutDir  = None
    print( "[timePsiMap] {0} is under processed...".format( args["job"] ), end="\n" )
    Output  = timePsiMap( job=args["job"], OutDir=OutDir, ksteps=args["ksteps"] )
    print( "\t [Done]" )
