import numpy               as np
import myStyle.cMap2D      as clm
import FetchEQU            as feq
import myStyle.LoadConfig  as lcf
import fLIB.fLIB__div2D    as div

def divBcheck( job     =None , InpDir  =None , CnsFile =None, \
               GridMode=None , x1Range =None , x2Range =None, \
               OutDir  =None , \
               silent  =False, config  =None  ):
    # --- [1] 引数チェック --- #
    if ( config   is None ): config   = lcf.LoadConfig()
    if ( job      is None ): job      = config["equ_job"]
    if ( GridMode is None ): GridMode = config["equ_GridMode"]
    if ( InpDir   is None ): InpDir   = "../job/{0}/{1}/"             .format( job, GridMode )
    if ( CnsFile  is None ): CnsFile  = "../job/{0}/{1}/parameter.dat".format( job, GridMode )
    if ( OutDir   is None ): OutDir   = "../job/{0}/{1}/pdf/"         .format( job, GridMode )
    # --- [2] データ呼び出し --- #
    Dkeys  = ["Bx","By","Bz"]
    Data   = feq.FetchEQU( job    =job    , InpDir  =InpDir  , CnsFile=CnsFile, \
                                keys   =Dkeys  , GridMode=GridMode, \
                                x1Range=x1Range, x2Range =x2Range , \
                                silent =silent , config  =config )
    # --- [3] プロット --- #
    xAxis  = Data["xAxis"]
    yAxis  = Data["yAxis"]
    divB   = div.div2d( x1g=Data["xAxis"], x2g=Data["yAxis"], \
                        u1 =Data["Bz"]   , u2 =Data["Bx"]   , \
                        coordinate="RTZ" , difftype="Backward_dx1")
    print( Data["Bx"].shape, Data["Bz"].shape, divB.shape )
    print( divB )

# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    print( "{0} is undir processing...".format( args["job"] ), end="" )
    cmp    = divBcheck( job =args["job"] )
    print( "\t [Done]" )
