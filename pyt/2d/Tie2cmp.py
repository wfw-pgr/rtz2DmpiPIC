import myStyle.cMap2D           as clm
import mPICsee.FetchPIC         as fpc
import myStyle.LoadConfig       as lcf
import myStyle.configSettings   as cfs
import strManipulate.keyInclude as stm
import Filter.nxLinearFilter2D  as flt

# ======================================== #
# ===  pic から 2次元カラーマップ描写  === #
# ======================================== #
def Tie2cmp( job  =None, jobDir=None, kstep=None, pngFile=None, pngDir=None, \
             keys =None, config=None ):
    # ------------------------------ #
    # --- [1] 引数チェック       --- #
    # ------------------------------ #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( pngDir  is None ): pngDir  = "{0}png/".format( jobDir )
    if ( pngFile is None ): pngFile = "{0}TiTe_{1}_.png".format( pngDir, job  )
    if ( keys    is None ): keys    = ["Te","Ti"]
    # ------------------------------ #
    # --- [2] データ読込         --- #
    # ------------------------------ #
    #  -- [2-1] key 設定         --  #
    config["pic_Temperature"] = True
    pkeys = [ key.replace( "T", "p" ) for key in keys ]
    if ( stm.keyInclude( pkeys, ["pi"] ) ):
        pkeys = pkeys + [ "pixx", "piyy", "pizz" ]
        pkeys.remove( "pi" ) 
    if ( stm.keyInclude( pkeys, ["pe"] ) ):
        pkeys = pkeys + [ "pexx", "peyy", "pezz" ]
        pkeys.remove( "pe" )
    if ( not( "psi" in keys )): pkeys   = pkeys + ["psi"]
    #  -- [2-2]   読込           --  #
    pData = fpc.FetchPIC( job =job  , jobDir=jobDir, kstep=kstep, keys=pkeys, config=config )
    xAxis = pData["xAxis"]
    yAxis = pData["yAxis"]
    flux  = pData["psi"  ]
    if ( stm.keyInclude( keys, ["Ti"] ) ): pData["pi"] = ( pData["pixx"] + pData["piyy"] + pData["pizz"] ) / 3.0
    if ( stm.keyInclude( keys, ["Te"] ) ): pData["pe"] = ( pData["pexx"] + pData["peyy"] + pData["pezz"] ) / 3.0
    #  -- [2-3] フィルタリング   --  #
    config["cmp_LinearFilt"] = 0.5
    config["cmp_nFilter"]    = 30
    pData = flt.nxLinearFilter2D( Data  =pData, keys  =pkeys, \
                                  x1Axis=xAxis, x2Axis=yAxis, coordinate="rtz", config=config )
    TData = {}
    for key in keys: TData[ key ] = pData[ key.replace("T","p") ]
    
    # ------------------------------ #
    # --- [3] プロット           --- #
    # ------------------------------ #
    cfs.configSettings( config=config, configType="cMap2D_thesis" )
    cfs.configSettings( config=config, configType="NoAxis" )
    config["clb_sw"]         = False
    config["cmp_ColorTable"] = "hot"
    config["cmp_Nlevels"]    = 255
    config["cmp_AutoLevel"]  = False
    config["cmp_MaxMin"]     = [0.0,0.10]
    config["MinimalOut"]     = True
    for key in keys:
        cMap    = TData[ key ]
        FigName = pngFile.replace( "_.png", "_{0}_{1:08}.png".format( key, kstep )  )
        clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis, cMap=cMap, Cntr=flux, \
                    FigName=FigName, config=config )



# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    pngDir = "png/Tie2cmp/"
    for kstep in args["ksteps"]:
        print( " [Tei2cmp] ( job , kstep ) = ( {0} , {1:8} )".format( args["job"], kstep ), end=""  )
        p2c    = Tie2cmp( job  =args["job"], jobDir=args["jobDir"], pngDir=pngDir, \
                          kstep=kstep      , keys  =args["key"]                    )
        print( "\t [Done]" )
