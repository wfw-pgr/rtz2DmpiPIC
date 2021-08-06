import numpy                  as np
import mPICsee.FetchPIC       as fpc
import myStyle.LoadConfig     as lcf
import myUtils.LoadConst      as lcn
import fLIB.fLIB__grad2d      as grd
import fLIB.fLIB__solveXOpt   as sxo
import myBasicAlgs.robustInv  as riv
import myStyle.configSettings as cfs
import mPICsee.Fetch1D        as f1d
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

    const  = lcn.LoadConst( InpFile=CnsFile )
    keys2D = ["uix", "uiy", "uiz" , "uex" , "uey", "uez" , \
              "pexy", "peyz", "ne" , "pixy", "piyz", "ni" ]
    Data2D = fpc.FetchPIC( job    =job    , jobDir =jobDir , kstep  =kstep, \
                           keys   =keys2D , config =config , x1Range=x1Range, x2Range=x2Range )
    xAxis  = + Data2D["xAxis"]
    yAxis  = + Data2D["yAxis"]
    keys1D = ["Ey", "uix", "uiy", "uiz" , "uex" , "uey", "uez" , "Bx", "By" , "Bz" ]
    Data1D = f1d.Fetch1D ( job    =job    , jobDir =jobDir , kstep  =kstep, \
                           keys   =keys1D , config =config , pltAxis=1    , x1Range=x1Range, x2Range=x2Range )
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
    Ey       = - Data1D["Ey"]
    uexB     = + Data1D["uez"]*Data1D["Bx"] - Data1D["uex"]*Data1D["Bz"]
    uixB     = + Data1D["uiz"]*Data1D["Bx"] - Data1D["uix"]*Data1D["Bz"]
    # -- (u.grad) u -- #
    grduet   = grd.grad2d( Data=Data2D["uey"], xg=xAxis, yg=yAxis  )
    grduit   = grd.grad2d( Data=Data2D["uiy"], xg=xAxis, yg=yAxis  )
    yInv     = riv.robustInv( yAxis )
    zI,rInv  = np.meshgrid( xAxis, yInv, indexing='xy' )
    ugrdue   = Data2D["uez"]*grduet["dfdx"] + Data2D["uex"]*grduet["dfdy"] + Data2D["uey"]* Data2D["uex" ] * rInv
    ugrdui   = Data2D["uiz"]*grduit["dfdx"] + Data2D["uix"]*grduit["dfdy"] + Data2D["uiy"]* Data2D["uix" ] * rInv
    ugrdue   = + np.ravel( np.mean( ugrdue*const["rme"], axis=meanAxis ) )
    ugrdui   = - np.ravel( np.mean( ugrdui*const["rmi"], axis=meanAxis ) )
    # -- div pe, pi -- #
    neInv    = riv.robustInv( Data2D["ne"], Flag__positive=True )
    niInv    = riv.robustInv( Data2D["ni"], Flag__positive=True )
    divPe_rt = ( grd.grad2d( Data=Data2D["pexy"], xg=xAxis, yg=yAxis ) )["dfdy"]
    divPe_zt = ( grd.grad2d( Data=Data2D["peyz"], xg=xAxis, yg=yAxis ) )["dfdx"]
    divPi_rt = ( grd.grad2d( Data=Data2D["pixy"], xg=xAxis, yg=yAxis ) )["dfdy"]
    divPi_zt = ( grd.grad2d( Data=Data2D["piyz"], xg=xAxis, yg=yAxis ) )["dfdx"]
    divPe    = + np.ravel( np.mean( (divPe_rt+divPe_zt)*neInv, axis=meanAxis ) )
    divPi    = - np.ravel( np.mean( (divPi_rt+divPi_zt)*niInv, axis=meanAxis ) )

    # ------------------------- #
    # --- [5] プロット      --- #
    # ------------------------- #
    pData   = { "Ey":Ey, "uexB":uexB, "ugrdue":ugrdue, "divPe":divPe, "uixB":uixB, "ugrdui":ugrdui, "divPi":divPi }
    label   = { "Ey":"E_y", \
                "uexB":r"$u_e \times B$", "ugrdue":r"$(u_e \cdot \nabla)u_e$", "divPe":r"$ \nabla \cdot P_{e}$", \
                "uixB":r"$u_i \times B$", "ugrdui":r"$(u_i \cdot \nabla)u_i$", "divPi":r"$ \nabla \cdot P_{i}$"  }
    cfs.configSettings( config=config, configType="plot1D_lateral" )
    FigName = OutFile.replace( ".png", "_e_{0:08}.png".format( kstep )  )
    # -- electron -- #
    fig     = pl1.plot1D( FigName=FigName, config=config )
    fig.addPlot( xAxis=Data2D[pAxiskey], yAxis =pData["Ey"]    , label=label["Ey"] )
    fig.addPlot( xAxis=Data2D[pAxiskey], yAxis =pData["uexB"]  , label=label["uexB"] )
    fig.addPlot( xAxis=Data2D[pAxiskey], yAxis =pData["ugrdue"], label=label["ugrdue"] )
    fig.addPlot( xAxis=Data2D[pAxiskey], yAxis =pData["divPe"] , label=label["divPe"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    # --   ion    -- #
    FigName = OutFile.replace( ".png", "_i_{0:08}.png".format( kstep )  )
    fig     = pl1.plot1D( FigName=FigName, config=config )
    fig.addPlot( xAxis=Data2D[pAxiskey], yAxis =pData["Ey"]    , label=label["Ey"] )
    fig.addPlot( xAxis=Data2D[pAxiskey], yAxis =pData["uixB"]  , label=label["uixB"] )
    fig.addPlot( xAxis=Data2D[pAxiskey], yAxis =pData["ugrdui"], label=label["ugrdui"] )
    fig.addPlot( xAxis=Data2D[pAxiskey], yAxis =pData["divPi"] , label=label["divPi"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()


# ------------------------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    OutDir  = "./png/ohmsLaw1D/"
    x1Range = [-8.0,+8.0]
    x2Range = [ 3.7, 6.7]
    pltAxis = 1
    for kstep in args["ksteps"]:
        ohmsLaw1D( job    =args["job"], jobDir =args["jobDir"], \
                   kstep  =kstep      , OutDir =OutDir        , \
                   pltAxis=pltAxis    , x1Range=x1Range       , x2Range=x2Range )
