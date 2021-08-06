import numpy                     as np
import myStyle.sfc3D             as sfc
import PICsees.FetchPIC          as fpc
import myStyle.LoadConfig        as lcf
import myConvert.where2D        as wh2
import fLIB.fLIB__LinearFilter2D as flt
import fLIB.fLIB__div2D          as div
import fLIB.fLIB__pyBiCGSTAB     as bcg

def poissonPotential( InpFile=None, InpDir=None, OutFile=None , OutDir=None, \
                      CnsFile=None, config=None ):
    # --- [1] 引数チェック --- #
    if ( config  is     None ): config   = lcf.LoadConfig()
    if ( CnsFile is     None ): CnsFile  = "../job/{0}/dat/constants.dat".format( config["job"] )
    if ( InpFile is     None ): InpFile  = "Field2D{0:08}.bin".format( config["kstep"] )
    if ( InpFile is     None ): InpDir   =   "../job/{0}/bin/".format( config["job"  ] )
    if ( OutFile is     None ): OutFile  = InpFile.replace( ".bin", ".png" )
    if ( OutDir  is not None ): OutFile  = OutDir + OutFile
    # --- [2] データ呼び出し --- #
    #  -- [2-1] Data Fetch  --  #
    keys        = ["Ex","Ey","Ez","psi"]
    Data        = fpc.FetchPIC( InpFile=InpFile, InpDir=InpDir, \
                                CnsFile=CnsFile, keys  =keys  , \
                                config =config )
    xAxis       = Data["xAxis"]
    yAxis       = Data["yAxis"]
    flux        = Data["psi"  ]
    #  -- [2-2] div E      -- #
    divE        = div.div2d( u1 =Data["Ez"]   , u2 =Data["Ex"]   , coordinate="RTZ", \
                             x1g=Data["xAxis"], x2g=Data["yAxis"], difftype  ="Central__dx2" )
    #   - divE 描きだし - #
    Flag__divEfigout= False
    if ( Flag__divEfigout ):
        FigName     = OutFile.replace( "Field2D", "divE" )
        fig         = sfc.sfc3D( xAxis  =yAxis,   yAxis =xAxis , sfc=divE, Cntr=flux, \
                                 FigName=FigName, config=config )
    #  -- 境界条件 -- #
    divE[ 0, :] = 0.0
    divE[-1, :] = 0.0
    divE[ :, 0] = 0.0
    divE[ :,-1] = 0.0
    rhs         = -divE
    phist       = bcg.pyBiCGSTAB( source=rhs, x1g=Data["xAxis"], x2g=Data["yAxis"], coordinate="RTZ", x2Min=0.0 )
    print( phist.shape, yAxis.shape, xAxis.shape )
    Data = wh2.where2D( Data=phist, x1Axis=yAxis, x2Axis=xAxis, \
                        x1Range=config["cmp_xRange"], x2Range=config["cmp_yRange"]  )
    flux = wh2.where2D( Data=flux, x1Axis=yAxis, x2Axis=xAxis, \
                        x1Range=config["cmp_xRange"], x2Range=config["cmp_yRange"]  )
    print( Data["Data"].shape, Data["x1Axis"].shape, Data["x2Axis"].shape )
    # --- [4] プロット  --- #
    config["cmp_AutoLevel"] = False
    config["cmp_MaxMin"]    = [-0.15,0.15]
    vmin, vmax              = -0.15, 0.15
    FigName                 = OutFile.replace( "Field2D", "sfc" )
    fig                     = sfc.sfc3D( xAxis  =Data["x1Axis"],   yAxis =Data["x2Axis"] ,
                                         sfc    =Data["Data"],  Cntr=flux["Data"], \
                                         elev   = 55 , azim=-80, \
                                         vmin   =vmin, vmax=vmax, \
                                         FigName=FigName, config=config )
    print( "[poissonPotential] {0} is outputed...".format( FigName ) )
    

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args      = rar.myRecvArgs()
    ConstFile = args["jobDir"] + "dat/constants.dat"
    InpDir    = args["jobDir"] + "bin/"
    OutDir    = args["jobDir"] + "png/"
    OutDir    = "./png/"
    config    = lcf.LoadConfig()

    # --- [3] コンフィグ --- #
    config["AutoRange"]       = False
    config["AutoTicks"]       = False
    config["FigSize"]         = (5,4)
    config["cmp_xRange"]      = [5.0,20.0]
    config["cmp_yRange"]      = [-5.0,+5.0]
    config["axes_x_ticks"]    = [5,10,15,20]
    config["axes_y_ticks"]    = [-5,0,+5]
    config["axes_x_off"]      = True
    config["axes_y_off"]      = True
    config["axes_z_off"]      = True
    config["cmp_Position"]    = [0.16,0.18,+0.80,0.89]
    config["cmp_AutoLevel"]   = True
    config["cmp_ColorTable"]  = "jet"
    config["clb_FontSize"]    = 12
    config["clb_Orientation"] = "vertical"
    config["clb_Position"]    = [0.86,0.20,0.90,0.90]
    config["clb_Title"]       = None
    config["MinimalOut"]      = True
    config["cmp_MaxMin"]      = [-0.003,0.003]

    
    for kstep in args["ksteps"]:
        binFile = "Field2D" + "{0:08d}".format( kstep ) + ".bin"
        vtk     = poissonPotential( InpFile=binFile  , InpDir=InpDir, OutDir=OutDir,
                                    CnsFile=ConstFile, config=config )
        print( binFile+" is processed..." )
