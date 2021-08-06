import numpy                     as np
import myStyle.cMap2D            as clm
import PICsees.FetchPIC          as fpc
import myStyle.LoadConfig        as lcf
import fLIB.fLIB__LinearFilter2D as flt
import fLIB.fLIB__gradInv2D      as gri

def staticPotential( InpFile=None, InpDir=None, OutFile=None , OutDir=None, \
                     CnsFile=None, config=None ):
    # --- [1] 引数チェック --- #
    if ( config  is     None ): config   = lcf.LoadConfig()
    if ( CnsFile is     None ): CnsFile  = "../job/{0}/dat/constants.dat".format( config["job"] )
    if ( InpFile is     None ): InpFile  = "Field2D{0:08}.bin".format( config["kstep"] )
    if ( InpFile is     None ): InpDir   =   "../job/{0}/bin/".format( config["job"  ] )
    if ( OutFile is     None ): OutFile  = InpFile.replace( ".bin", ".png" )
    if ( OutDir  is not None ): OutFile  = OutDir + OutFile
    # --- [2] データ呼び出し --- #
    keys        = ["Ex","Ey","Ez","psi"]
    Data        = fpc.FetchPIC( InpFile=InpFile, InpDir=InpDir, \
                                CnsFile=CnsFile, keys  =keys  , \
                                config =config )
    xAxis       = Data["xAxis"]
    yAxis       = Data["yAxis"]
    flux        = Data["psi"  ]
    #  -- [2-2] gradInv      -- #
    phist       = gri.gradInv2D( E1 = - Data["Ez"]   , E2 = - Data["Ex"], \
                                 x1g=   Data["xAxis"], x2g=   Data["yAxis"], \
                                 startPoint="center" )
    print( xAxis.shape, yAxis.shape, phist.shape )
    
    # --- [4] プロット  --- #
    FigName                 = OutFile.replace( "Field2D", "staticP" )
    fig                     = clm.cMap2D( xAxis  =yAxis,   yAxis =xAxis , cMap=phist, Cntr=flux, \
                                          FigName=FigName, config=config )
    print( "[staticPotential] {0} is outputed...".format( FigName ) )
    

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
    config["FigSize"]         = (3.0,5.0)
    config["AutoRange"]       = True
    config["cmp_xRange"]      = [5.0,20.0]
    config["cmp_yRange"]      = [-5.0,+5.0]
    config["AutoTicks"]       = True
    config["axes_x_ticks"]    = [5,10,15,20]
    config["axes_y_ticks"]    = [-5,0,+5]
    config["cmp_Position"]    = [0.16,0.18,+0.80,0.89]
    config["cmp_AutoLevel"]   = True
    config["cmp_ColorTable"]  = "plasma"
    config["clb_orientation"] = "vertical"
    config["clb_title"]       = None
    config["clb_title_pos"]   = [0.85,0.91]
    config["clb_title_size"]  = 24
    config["clb_position"]    = [0.92,0.20,0.96,0.89]
    config["MinimalOut"]      = False
    config["cmp_MaxMin"]      = [-0.003,0.003]

    
    for kstep in args["ksteps"]:
        binFile = "Field2D" + "{0:08d}".format( kstep ) + ".bin"
        vtk     = staticPotential( InpFile=binFile, InpDir=InpDir, OutDir=OutDir, CnsFile=ConstFile, config=config )
        print( binFile+" is processed..." )
