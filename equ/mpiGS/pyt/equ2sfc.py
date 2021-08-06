import numpy               as np
import myStyle.sfc3D       as sfc
import FetchEQU            as feq
import myStyle.LoadConfig  as lcf

def equ2sfc( job     =None , \
             InpDir  =None , CnsFile =None, \
             keys    =None , GridMode=None, \
             x1Range =None , x2Range =None, \
             OutDir  =None , \
             silent  =False, config  =None    ):
    # --- [1] 引数チェック --- #
    if ( config   is None ): config   = lcf.LoadConfig()
    if ( job      is None ): job      = config["equ_job"]
    if ( GridMode is None ): GridMode = config["equ_GridMode"]
    if ( keys     is None ): keys     = [ config["equ_key"], "psi" ]
    if ( InpDir   is None ): InpDir   = "../job/{0}/{1}/"             .format( job, GridMode )
    if ( CnsFile  is None ): CnsFile  = "../job/{0}/{1}/parameter.dat".format( job, GridMode )
    if ( OutDir   is None ): OutDir   = "../job/{0}/{1}/pdf/"         .format( job, GridMode )
    # --- [2] データ呼び出し --- #
    Dkeys       = keys
    if ( not( "psi" in keys ) ): Dkeys =keys + ["psi"]
    Data        = feq.FetchEQU( job    =job    , InpDir  =InpDir  , CnsFile=CnsFile, \
                                keys   =Dkeys  , GridMode=GridMode, \
                                x1Range=x1Range, x2Range =x2Range , \
                                silent =silent , config  =config )
    # --- [3] プロット --- #
    xAxis       = + Data["xAxis"]
    yAxis       = + Data["yAxis"]
    flux        = + Data["psi"  ]
    for key in keys:
        sfcMap  = Data[ key ]
        FigName = "{0}sfc-{1}_.png".format( OutDir, key )
        print( "[equ2sfc] {0} is under processing...".format( FigName ), end="" )
        fig     = sfc.sfc3D( xAxis  =yAxis  , yAxis =xAxis , sfc=sfcMap, Cntr=flux, \
                             FigName=FigName, config=config )
        print( "\t [Done]" )


# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    OutDir = "./png/"
    print( "{0} is undir processing...".format( args["job"] ) )
    keys   = [ "Bx" , "By" , "Bz" , "Jx", "Jy", "Jz",
               "psi", "rho", "prs", "ux", "uy", "uz"  ]
    sfc    = equ2sfc( job =args["job"], OutDir=OutDir, keys=keys )
