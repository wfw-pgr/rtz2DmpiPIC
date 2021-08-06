import numpy                   as np
import myStyle.LoadConfig      as lcf
import myUtils.LoadConst       as lcn
import mPICsee.FetchPIC        as fpc
import myStyle.cMap2D          as clm
import myStyle.configSettings  as cfs
import myAnalysis.createMask   as msk
import Filter.nxLinearFilter2D as flt
import myBasicAlgs.cylindricalVolume as cvl

# ============================================ #
# ===  イオン / 電子温度 - 質量比 依存性   === #
# ============================================ #
def mr_Upt( job=None, kstep=None, jobDir=None, datDir=None, config=None, pngDir=None, Flag__Display=False, CnsFile=None ):
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}".format( config["pic_jobDir"], job )
    if ( datDir  is None ): datDir  = "{0}dat/".format( jobDir )
    if ( pngDir  is None ): pngDir  = "{0}png/".format( jobDir )
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )

    # ------------------------------ #
    # --- [2] データ読込         --- #
    # ------------------------------ #
    config["pic_Temperature"] = True
    const      = lcn.LoadConst( InpFile=CnsFile )
    Pkeys      = ["pexx","peyy","pezz","pixx","piyy","pizz","ne","ni","psi",
                  "uix","uiy","uiz","uex","uey","uez" ]
    Data       = fpc.FetchPIC( job =job  , jobDir=jobDir, kstep =kstep, \
                               keys=Pkeys, config=config, RowData=True )
    xAxis      =   Data["xAxis"]
    yAxis      =   Data["yAxis"]
    flux       =   Data["psi"  ]
    ethi       = 0.5 *              ( Data["pixx"]   + Data["piyy"]   + Data["pizz"]   )
    ethe       = 0.5 *              ( Data["pexx"]   + Data["peyy"]   + Data["pezz"]   )
    eFli       = 0.5 * Data["ni"] * ( Data["uix"]**2 + Data["uiy"]**2 + Data["uiz"]**2 ) 
    eFle       = 0.5 * Data["ne"] * ( Data["uex"]**2 + Data["uey"]**2 + Data["uez"]**2 ) 
    mask       = msk.createMask( Flux=Data["psi"], Flag__LCFSMode=True )
    volg       = cvl.cylindricalVolume( zAxis=xAxis, rAxis=yAxis )
    volmSum    = np.sum(               volg * mask )
    etheAvg    = np.sum(  ethe       * volg * mask ) / volmSum
    ethiAvg    = np.sum(  ethi       * volg * mask ) / volmSum
    eFleAvg    = np.sum(  eFle       * volg * mask ) / volmSum
    eFliAvg    = np.sum(  eFli       * volg * mask ) / volmSum
    epteAvg    = np.sum( (ethe+eFle) * volg * mask ) / volmSum
    eptiAvg    = np.sum( (ethi+eFli) * volg * mask ) / volmSum
    print( "[mr_Upt] results ----" )
    print( " volm : {0}".format( volmSum ) )
    print( " ethe : {0}".format( etheAvg ) )
    print( " ethi : {0}".format( ethiAvg ) )
    print( " eFle : {0}".format( eFleAvg ) )
    print( " eFli : {0}".format( eFliAvg ) )
    print( " epte : {0}".format( epteAvg ) )
    print( " epti : {0}".format( eptiAvg ) )
    return( [ float(kstep), const["mr"], volmSum, etheAvg, ethiAvg, eFleAvg, eFliAvg, epteAvg, eptiAvg ] )

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
    ret     = mr_Upt( job=args["job"], jobDir=args["jobDir"], kstep=args["kstep"], pngDir=pngDir )
    print( ret )
