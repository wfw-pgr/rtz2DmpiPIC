import numpy                  as np
import mPICsee.FetchPIC       as fpc
import myStyle.LoadConfig     as lcf
import myUtils.LoadConst      as lcn
import fLIB.fLIB__grad2d      as grd
import fLIB.fLIB__solveXOpt   as sxo
import myBasicAlgs.robustInv  as riv
import myStyle.configSettings as cfs
import myBasicAlgs.robustInv  as inv
import myBasicAlgs.makegrid   as mkg
import myStyle.plot1D         as pl1
    
def ohmsLaw1D( job  =None, jobDir=None, kstep  =None, OutFile=None, OutDir =None, CnsFile=None, \
               keys =None, config=None, pltAxis=None, x1Range=None, x2Range=None, Flag__showXptPosition=True ):
    # ------------------------ #
    # --- [1] 引数チェック --- #
    # ------------------------ #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}png/".format( jobDir )
    if ( OutFile is None ): OutFile = "OhmsLaw1D.png"
    if ( pltAxis is None ): pltAxis = config["plt_Axis"]
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    OutFile = OutDir + OutFile
    
    # -------------------------- #
    # --- [2] データ呼び出し --- #
    # -------------------------- #
    #  --  軸設定  -- #
    if   ( pltAxis == 1 ): # -- z-Axis -- #
        pAxiskey = "xAxis"
        meanAxis = 0
        if ( x1Range is None ): x1Range = config["Data_x1Range" ]
        if ( x2Range is None ): x2Range = config["Data_avgRange"]
    elif ( pltAxis == 2 ): # -- r-Axis -- #
        pAxiskey = "yAxis"
        meanAxis = 1
        if ( x1Range is None ): x1Range = config["Data_avgRange"]
        if ( x2Range is None ): x2Range = config["Data_x2Range" ]
    #  --  データ  -- #
    const   = lcn.LoadConst( InpFile=CnsFile )
    LFalpha = 0.1
    nFilter = 10
    keys    = ["Ey"  , "Bx"  , "By", "Bz"  , "uex" , "uey", "uez" , "uix", "uiy", "uiz" , \
               "pexy", "peyz", "ne", "pixy", "piyz", "ni" ]
    Data    = fpc.FetchPIC( job       =job , jobDir    =jobDir , kstep  =kstep, \
                            keys      =keys, config    =config , x1Range=x1Range, x2Range=x2Range, \
                            FilterKeys=keys, LinearFilt=LFalpha, nFilter=nFilter, RowData=True     )
    #  --  座標軸  -- #
    xAxis   = Data["xAxis"]
    yAxis   = Data["yAxis"]
    yInv    = inv.robustInv( yAxis )
    yIg     = ( mkg.makegrid( x1=yInv, x2=xAxis, indexing='xy' ) )["x1g"]
    qe, qi  = -1.0, +1.0
    #  --  X点位置 -- #
    if ( Flag__showXptPosition ):
        config["phys_LengUnit"] = "c/wce"
        Flux    = fpc.FetchPIC( job=job, jobDir=jobDir, kstep=kstep, keys=["psi"], config =config )
        ptOXO   = ( sxo.solveXOpt( psi=Flux["psi"]) )["ptOXO"]
        x1Xpt, x2Xpt = ( Flux["xAxis"] )[ ptOXO[0,1] ], ( Flux["yAxis"] )[ ptOXO[1,1] ]
        print( "[OhmsLaw1D] Xpt position (x1,x2) = ( {0},{1} )".format( x1Xpt, x2Xpt ) )
    
    # ------------------------------ #
    # --- [3] Ohm's Law ( 電子 ) --- #
    # ------------------------------ #
    # --  Ey, uxB   -- #
    Ey       = - np.ravel( np.mean( Data["Ey"], axis=meanAxis ) )
    uexB     =   Data["uez"]*Data["Bx"] - Data["uex"]*Data["Bz"]
    uixB     =   Data["uiz"]*Data["Bx"] - Data["uix"]*Data["Bz"]
    uexB     =   np.ravel( np.mean( uexB      , axis=meanAxis ) )
    uixB     =   np.ravel( np.mean( uixB      , axis=meanAxis ) )
    # -- (u.grad) u -- #
    grduet   = grd.grad2d( Data=Data["uey"], xg=xAxis, yg=yAxis  )
    grduit   = grd.grad2d( Data=Data["uiy"], xg=xAxis, yg=yAxis  )
    ugrdue   = Data["uez"]*grduet["dfdx"] + Data["uex"]*grduet["dfdy"] + Data["uey"]*Data["uex" ]*yIg
    ugrdui   = Data["uiz"]*grduit["dfdx"] + Data["uix"]*grduit["dfdy"] + Data["uiy"]*Data["uix" ]*yIg
    ugrdue   = - np.ravel( np.mean( ugrdue, axis=meanAxis ) ) * ( const["rme"] / qe )
    ugrdui   = - np.ravel( np.mean( ugrdui, axis=meanAxis ) ) * ( const["rme"] / qi )
    # -- div pe, pi -- #
    qneInv    = riv.robustInv( Data["ne"], Flag__positive=True ) / qe
    qniInv    = riv.robustInv( Data["ni"], Flag__positive=True ) / qi
    divPe_rt = ( grd.grad2d( Data=Data["pexy"], xg=xAxis, yg=yAxis ) )["dfdy"]
    divPe_zt = ( grd.grad2d( Data=Data["peyz"], xg=xAxis, yg=yAxis ) )["dfdx"]
    divPi_rt = ( grd.grad2d( Data=Data["pixy"], xg=xAxis, yg=yAxis ) )["dfdy"]
    divPi_zt = ( grd.grad2d( Data=Data["piyz"], xg=xAxis, yg=yAxis ) )["dfdx"]
    divPe    = - np.ravel( np.mean( (divPe_rt+divPe_zt)*qneInv, axis=meanAxis ) )
    divPi    = - np.ravel( np.mean( (divPi_rt+divPi_zt)*qniInv, axis=meanAxis ) )

    # ------------------------- #
    # --- [5] プロット      --- #
    # ------------------------- #
    pData   = { "Ey":Ey, "uexB":uexB, "ugrdue":ugrdue, "divPe":divPe, "uixB":uixB, "ugrdui":ugrdui, "divPi":divPi }
    label   = { "Ey":"$E_y$", \
                "uexB":r"$u_e \times B$", "ugrdue":r"$(u_e \cdot \nabla)u_e$", "divPe":r"$ \nabla \cdot P_{e}$", \
                "uixB":r"$u_i \times B$", "ugrdui":r"$(u_i \cdot \nabla)u_i$", "divPi":r"$ \nabla \cdot P_{i}$"  }
    cfs.configSettings( config=config, configType="plot1D_lateral" )
    config["plt_MedianFilt"] = 0
    config["plt_yAutoRange"] = False
    config["plt_yRange"]     = [-0.1,0.1]
    config["cursor_y"]       = [0.0]
    config["plt_GaussFilt"]  = 2.0
    config["plt_LinearFilt"] = 0.0
    config["plt_LegNColumn"] = 2
    FigName = OutFile.replace( ".png", "_e_{0:08}.png".format( kstep )  )
    # -- electron -- #
    fig     = pl1.plot1D( FigName=FigName, config=config )
    fig.addPlot( xAxis=Data[pAxiskey], yAxis =pData["Ey"]    , label=label["Ey"]     )
    fig.addPlot( xAxis=Data[pAxiskey], yAxis =pData["uexB"]  , label=label["uexB"]   )
    fig.addPlot( xAxis=Data[pAxiskey], yAxis =pData["ugrdue"], label=label["ugrdue"] )
    fig.addPlot( xAxis=Data[pAxiskey], yAxis =pData["divPe"] , label=label["divPe"]  )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    # --   ion    -- #
    FigName = OutFile.replace( ".png", "_i_{0:08}.png".format( kstep )  )
    fig     = pl1.plot1D( FigName=FigName, config=config )
    fig.addPlot( xAxis=Data[pAxiskey], yAxis =pData["Ey"]    , label=label["Ey"]     )
    fig.addPlot( xAxis=Data[pAxiskey], yAxis =pData["uixB"]  , label=label["uixB"]   )
    fig.addPlot( xAxis=Data[pAxiskey], yAxis =pData["ugrdui"], label=label["ugrdui"] )
    fig.addPlot( xAxis=Data[pAxiskey], yAxis =pData["divPi"] , label=label["divPi"]  )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()


# ------------------------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    OutDir  = "./png/ohmsLaw1D/"
    x1Range = [-8.0,+8.0]
    x2Range = [ 4.5, 4.7]
    pltAxis = 1
    for kstep in args["ksteps"]:
        ohmsLaw1D( job    =args["job"], jobDir =args["jobDir"], \
                   kstep  =kstep      , OutDir =OutDir        , \
                   pltAxis=pltAxis    , x1Range=x1Range       , x2Range=x2Range )
