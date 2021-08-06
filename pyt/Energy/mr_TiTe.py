import numpy                   as np
import myStyle.LoadConfig      as lcf
import mPICsee.FetchPIC        as fpc
import myStyle.cMap2D          as clm
import myStyle.configSettings  as cfs
import myAnalysis.createMask   as msk
import Filter.nxLinearFilter2D as flt

# ============================================ #
# ===  イオン / 電子温度 - 質量比 依存性   === #
# ============================================ #
def mr_TiTe( job=None, kstep=None, jobDir=None, datDir=None, config=None, pngDir=None, Flag__Display=False ):
    if ( config is None ): config = lcf.LoadConfig()
    if ( job    is None ): job    = config["pic_job"]
    if ( kstep  is None ): kstep  = config["pic_kstep"]
    if ( jobDir is None ): jobDir = "{0}{1}".format( config["pic_jobDir"], job )
    if ( datDir is None ): datDir = "{0}dat/".format( jobDir )
    if ( pngDir is None ): pngDir = "{0}png/".format( jobDir )

    # ------------------------------ #
    # --- [2] データ読込         --- #
    # ------------------------------ #
    config["pic_Temperature"] = True
    Pkeys      = [ "pexx", "peyy", "pezz", "pixx", "piyy", "pizz", "ne", "ni", "psi" ]
    Data       = fpc.FetchPIC( job =job  , jobDir=jobDir, kstep =kstep, \
                               keys=Pkeys, config=config, RowData=True )
    xAxis      =   Data["xAxis"]
    yAxis      =   Data["yAxis"]
    flux       =   Data["psi"  ]
    Data["Ti"] = ( Data["pixx"] + Data["piyy"] + Data["pizz"] ) / 3.0
    Data["Te"] = ( Data["pexx"] + Data["peyy"] + Data["pezz"] ) / 3.0
    mask       = msk.createMask( Flux=Data["psi"], Flag__LCFSMode=True )
    import myBasicAlgs.cylindricalVolume as cvl
    volg       = cvl.cylindricalVolume( zAxis=xAxis, rAxis=yAxis )
    volmSum    = np.sum(              volg * mask )
    etheAvg    = np.sum( Data["Te"] * volg * mask ) / volmSum
    ethiAvg    = np.sum( Data["Ti"] * volg * mask ) / volmSum
    print( volg.shape )
    print( " volm : {0}".format( volmSum ) )
    print( " ethe : {0}".format( etheAvg ) )
    print( " ethi : {0}".format( ethiAvg ) )

    # ---------------------------------------- #
    # --- [3] ディスプレイ                 --- #
    # ---------------------------------------- #
    if ( Flag__Display ):
        config["cmp_LinearFilt"] = 0.2
        config["cmp_nFilter"]    = 10
        ret        = flt.nxLinearFilter2D( Data  =Data , keys=["Ti","Te"], \
                                           x1Axis=xAxis, x2Axis=yAxis, coordinate="rtz", config=config )
        Data["Ti"], Data["Te"] = ret["Ti"], ret["Te"]
        cfs.configSettings( config=config, configType="cMap2D_thesis" )
        for key in ["Ti","Te"]:
            FigName = "./{0}Masked_{1:08}.png".format( key, kstep )
            clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis, cMap=Data[key], Cntr=flux, \
                        FigName=FigName, config=config )


# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    pngDir  = "./"
    ret     = mr_TiTe( job=args["job"], jobDir=args["jobDir"], kstep=args["kstep"], pngDir=pngDir )
