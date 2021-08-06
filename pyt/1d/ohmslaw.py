import numpy               as np
import myStyle.plot1D      as pl1
import PICsees.FetchPIC    as fpc
import myStyle.LoadConfig  as lcf
import fLIB.fLIB__grad2d   as grd
import myUtils.LoadConst   as lcn

def ohmslaw( InpFile=None, InpDir =None, \
             OutFile=None, OutDir =None, \
             CnsFile=None, config =None, \
             labels =None, pltAxis=None, \
             x1Range=None, x2Range=None  ):
    # --- [1] 引数チェック --- #
    if ( config  is     None ): config  = lcf.LoadConfig()
    if ( CnsFile is     None ): CnsFile = "../job/{0}/dat/constants.dat".format( config["job"] )
    if ( InpFile is     None ): InpFile = "Field2D{0:08}.bin".format( config["kstep"] )
    if ( InpDir  is     None ): InpDir  =   "../job/{0}/bin/".format( config["job"  ] )
    if ( OutFile is     None ): OutFile = InpFile.replace( ".bin", ".png" )
    if ( OutDir  is not None ): OutFile = OutDir + OutFile
    if ( pltAxis is     None ): pltAxis = config["plt_Axis"]
    
    # --- [2] データ呼び出し --- #
    #  -- 軸の選択 / 平均化  --  #
    if   ( pltAxis == 1 ):
        meanAxis = 0
        if ( x1Range is None ): x1Range = config["Data_x1Range" ]
        if ( x2Range is None ): x2Range = config["Data_avgRange"]
    elif ( pltAxis == 2 ):
        meanAxis = 1
        if ( x1Range is None ): x1Range = config["Data_avgRange"]
        if ( x2Range is None ): x2Range = config["Data_x2Range" ]
    else:
        sys.exit( "[ERROR] pst_Axis?? [ERROR]" )
    #  -- データ読み出し -- #
    keys   = ["Ey","uix","uiy","uiz","uex","uey","uez", "psi",
              "Bx","By","Bz","pexy","peyz", "ne", "pixy","piyz", "ni" ]
    Field  = fpc.FetchPIC( InpFile=InpFile, InpDir=InpDir, \
                           CnsFile=CnsFile, keys  =keys  , \
                           config =config , \
                           x1Range=x1Range, x2Range=x2Range )
    const  = lcn.LoadConst( InpFile=CnsFile )
    # ---   ion    --- #
    ion__Analysis = True
    if ( ion__Analysis ):
        # -- Ey  -- #
        Ey     = - np.ravel( np.mean( Field["Ey"]                                        , axis=meanAxis ) )
        # -- uxB -- #
        uixB   = + np.ravel( np.mean( Field["uiz"]*Field["Bx"] - Field["uix"]*Field["Bz"], axis=meanAxis ) )
        # -- (ui.grad) ui -- #
        grduit = grd.grad2d ( Data=Field["uiy"], xg=Field["xAxis"], yg=Field["yAxis"] )
        zI,rInv= np.meshgrid( Field["xAxis"], 1.0 / ( Field["yAxis"] ), indexing='xy' )
        ugrdu  = Field["uiz"]*grduit["dfdx"] + Field["uix"]*grduit["dfdy"] + Field["uiy"]* Field["uix" ] * rInv
        ugrdu  = - np.ravel( np.mean( ugrdu*const["rmi"], axis=meanAxis ) )
        print( const["rmi"] )
        # -- div pi -- #
        niInv  = 1.0 / ( Field["ni"] )
        divPrt = ( grd.grad2d( Data=Field["pixy"], xg=Field["xAxis"], yg=Field["yAxis"] ) )["dfdy"]
        divPzt = ( grd.grad2d( Data=Field["piyz"], xg=Field["xAxis"], yg=Field["yAxis"] ) )["dfdx"]
        divPi  = - np.ravel( np.mean( (divPrt+divPzt)*niInv, axis=meanAxis ) )
        # -- Residual -- #
        res    = - Ey + uixB + ugrdu + divPi

        # --- [3] プロット --- #
        OutFile= OutDir + "iOhmsLaw.png"
        labels = [r"$-E_y$", r"$u_i \times B$" , r"$m_i / e (u_i \cdot \nabla) u_i$",
                  r"$\nabla \cdot P_i / e n_i$", "Residual"]
        fig    = pl1.plot1D( FigName=OutFile, config=config )
        xAxis  = Field["xAxis"] if ( pltAxis == 1 ) else Field["yAxis"]
        fig.addPlot( xAxis=xAxis, yAxis=Ey,    label=labels[0] )
        fig.addPlot( xAxis=xAxis, yAxis=uixB,  label=labels[1] )
        fig.addPlot( xAxis=xAxis, yAxis=ugrdu, label=labels[2] )
        fig.addPlot( xAxis=xAxis, yAxis=divPi, label=labels[3] )
        # fig.addPlot( xAxis=xAxis, yAxis=res,   label=labels[4] )
        fig.addLegend()
        fig.setAxis()
        fig.writeFigure()
        print( "[ohmslaw] {0} is outputed...".format( OutFile ) )

        
    # --- Electron --- #
    Electron__Analysis = True
    if ( Electron__Analysis ):
        # -- Ey  -- #
        Ey     = - np.ravel( np.mean( Field["Ey"]                                        , axis=meanAxis ) )
        # -- uxB -- #
        uexB   = + np.ravel( np.mean( Field["uez"]*Field["Bx"] - Field["uex"]*Field["Bz"], axis=meanAxis ) )
        # -- (ue.grad) ue -- #
        grduet = grd.grad2d( Data=Field["uey"], xg=Field["xAxis"], yg=Field["yAxis"]  )
        zI,rInv= np.meshgrid( Field["xAxis"], 1.0 / ( Field["yAxis"] ), indexing='xy' )
        ugrdu  = Field["uez"]*grduet["dfdx"] + Field["uex"]*grduet["dfdy"] + Field["uey"]* Field["uex" ] * rInv
        ugrdu  = + np.ravel( np.mean( ugrdu, axis=meanAxis ) )
        # -- div pe -- #
        neInv  = 1.0 / ( Field["ne"] )
        divPrt = ( grd.grad2d( Data=Field["pexy"], xg=Field["xAxis"], yg=Field["yAxis"] ) )["dfdy"]
        divPzt = ( grd.grad2d( Data=Field["peyz"], xg=Field["xAxis"], yg=Field["yAxis"] ) )["dfdx"]
        divPe  = + np.ravel( np.mean( (divPrt+divPzt)*neInv, axis=meanAxis ) )
        # -- Residual -- #
        res    = - Ey + uexB + ugrdu + divPe
        
        # --- [3] プロット --- #
        OutFile= OutDir + "eOhmsLaw.png"
        labels = [r"$-E_y$", r"$u_e \times B$" , r"$m_e / e (u_e \cdot \nabla) u_e$",
                  r"$\nabla \cdot P_e / e n_e$", "Residual"]
        fig    = pl1.plot1D( FigName=OutFile, config=config )
        xAxis  = Field["xAxis"] if ( pltAxis == 1 ) else Field["yAxis"]
        fig.addPlot( xAxis=xAxis, yAxis=Ey,    label=labels[0] )
        fig.addPlot( xAxis=xAxis, yAxis=uexB,  label=labels[1] )
        fig.addPlot( xAxis=xAxis, yAxis=ugrdu, label=labels[2] )
        fig.addPlot( xAxis=xAxis, yAxis=divPe, label=labels[3] )
        fig.addPlot( xAxis=xAxis, yAxis=res,   label=labels[4] )
        fig.addLegend()
        fig.setAxis()
        fig.writeFigure()
        print( "[ohmslaw] {0} is outputed...".format( OutFile ) )



# --------------------------------------- #
if ( __name__=="__main__" ):
    # --- 引数の設定 --- #
    import myUtils.myRecvArgs as rar
    args      = rar.myRecvArgs()
    CnsFile   = args["jobDir"] + "dat/constants.dat"
    InpDir    = args["jobDir"] + "bin/"
    OutDir    = "./png/"
    config    = lcf.LoadConfig()

    #  -- 準備 psiX の位置を知る -- #
    inquiry_PsiX = True
    if ( inquiry_PsiX ):
        kstep   = args["ksteps"][0]
        InpFile = "Field2D{0:08}.bin".format( kstep )
        x1Range = [-3.,+3.]
        x2Range = [+6.,18.]
        import PICsees.getXpt_FromPIC as gPX
        ret = gPX.getXpt_FromPIC( InpFile=InpFile, InpDir =InpDir, \
                                  CnsFile=CnsFile, config =config, \
                                  x1Range=x1Range, x2Range=x2Range \
                              )
        print()
        print( "psiX    ::  {0}     ".format( ret[0] ) )
        print( "(iX,jX) :: ({0},{1})".format( ret[1],ret[2] ) )
        print( "(xX,yX) :: ({0},{1})".format( ret[3],ret[4] ) )
        print()
    #  -- psi を描いて 知る --  #
    check_PsiX = True
    if ( check_PsiX ):
        kstep   = args["ksteps"][0]
        InpFile = "Field2D{0:08}.bin".format( kstep )
        x1Range = [-2.,+2.]
        x2Range = [6.0,18.0]
        Field   = fpc.FetchPIC( InpFile=InpFile, InpDir =InpDir , CnsFile=CnsFile, config =config, \
                                keys   =["psi"], x1Range=x1Range, x2Range=x2Range )
        import myStyle.cMap2D as clm
        fig     = clm.cMap2D( xAxis=Field["yAxis"], yAxis=Field["xAxis"], cMap=Field["psi"] )
        
    
    # --- コンフィグの設定 --- #
    config["FigSize"]         = (5,3)
    config["AutoRange"]       = False
    config["AutoTicks"]       = True
    config["Data_x1Range"]    = [-10,+10]
    config["Data_avgRange"]   = [+12.5,13.5]
    config["plt_xRange"]      = [-5.,+5.]
    config["plt_yRange"]      = [-0.05,+0.05]
    config["plt_Axis"]        = 1
    config["plt_Position"]    = [0.22,0.22,0.95,0.95]
    config["plt_LegNColumn"]  = 2
    config["plt_LegFontSize"] = 10
    config["plt_LegLocation"] = "lower-left"
    config["axes_x_Nticks"]   = 5
    config["axes_y_Nticks"]   = 5
    config["xTitle"]          = "Z"
    config["yTitle"]          = ""
    config["MinimalOut"]      = False
    # --- 描きだし --- #
    for kstep in args["ksteps"]:
        binFile = "Field2D{0:08d}.bin".format( kstep )
        vtk     = ohmslaw( InpFile=binFile, InpDir=InpDir, OutDir=OutDir, CnsFile=CnsFile, config=config )
        print( binFile+" is processed..." )
