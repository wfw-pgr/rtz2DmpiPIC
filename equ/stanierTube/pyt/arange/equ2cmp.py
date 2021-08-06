import myStyle.cMap2D         as clm
import mPICsee.FetchEQU       as feq
import myStyle.configSettings as cfs
import myStyle.LoadConfig     as lcf
import myUtils.LoadConst      as lcn

def equ2cmp( job     =None , InpDir  =None , CnsFile =None, keys    =None , GridMode=None, \
             x1Range =None , x2Range =None , OutDir  =None, Flag__electronScale=True, \
             silent  =False, config  =None ):
    # --- [1] 引数チェック --- #
    if ( config   is None ): config   = lcf.LoadConfig()
    if ( job      is None ): job      = config["equ_job"]
    if ( GridMode is None ): GridMode = config["equ_GridMode"]
    if ( keys     is None ): keys     = [ config["equ_key"], "psi" ]
    if ( InpDir   is None ): InpDir   = "../../job/{0}/{1}/".format( job, GridMode )
    if ( CnsFile  is None ): CnsFile  = "{0}parameter.dat"  .format( InpDir )
    if ( OutDir   is None ): OutDir   = "{0}pdf/"           .format( InpDir )
    # --- [2] データ呼び出し --- #
    const       = lcn.LoadConst( InpFile=CnsFile )
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
    if ( Flag__electronScale ):
        xAxis = xAxis * const["wpewce"]
        yAxis = yAxis * const["wpewce"]
    cfs.configSettings( config=config, configType="cMap2D_thesis" )
    config["axes_x_NoLabel"] = True
    config["axes_y_NoLabel"] = True
    # --  Jy -- #
    key = "Jy"
    config["cnt_Color"]      = "grey"
    config["cmp_ColorTable"] = "jet"
    config["cmp_AutoLevel"]  = False
    config["cmp_MaxMin"]     = [0,0.2]
    FigName = "{0}{1}_.png".format( OutDir, key )
    print( "[equ2cmp] {0:30} is under processing...".format( FigName ), end="" )
    clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis , cMap=Data[key], Cntr=flux, \
                FigName=FigName, config=config )
    print( "\t [Done]" )
    # --  By -- #
    key = "By"
    config["cnt_Color"]      = "grey"
    config["cmp_ColorTable"] = "bwr"
    config["cmp_AutoLevel"]  = False
    config["cmp_MaxMin"]     = [-5,5]
    FigName = "{0}{1}_.png".format( OutDir, key )
    print( "[equ2cmp] {0:30} is under processing...".format( FigName ), end="" )
    clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis , cMap=Data[key], Cntr=flux, \
                FigName=FigName, config=config )
    print( "\t [Done]" )
    # -- psi -- #
    key = "psi"
    config["cnt_Color"]      = "DeepSkyBlue"
    config["cmp_ColorTable"] = "hot"
    config["cmp_AutoLevel"]  = False
    config["cmp_MaxMin"]     = [-1.4,1.0]
    FigName = "{0}{1}_.png".format( OutDir, key )
    print( "[equ2cmp] {0:30} is under processing...".format( FigName ), end="" )
    clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis , cMap=Data[key], Cntr=flux, \
                FigName=FigName, config=config )
    print( "\t [Done]" )
    

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    OutDir = "./png/"
    print( "{0} is undir processing...".format( args["job"] ) )
    keys   = [ "Jy", "psi", "By" ]
    cmp    = equ2cmp( job =args["job"], OutDir=OutDir, keys=keys )
