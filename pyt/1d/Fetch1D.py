import numpy              as np
import mPICsee.FetchPIC   as fpc
import myStyle.LoadConfig as lcf

def Fetch1D( job    =None, jobDir =None, kstep  =None, keys =None, config=None, \
             pltAxis=None, x1Range=None, x2Range=None ):
    # ------------------------ #
    # --- [1] 引数チェック --- #
    # ------------------------ #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( pltAxis is None ): pltAxis = config["plt_Axis"]

    # -------------------------- #
    # --- [2] データ呼び出し --- #
    # -------------------------- #
    if   ( pltAxis == 1 ):
        pAxiskey = "xAxis"
        meanAxis = 0
        if ( x1Range is None ): x1Range = config["Data_x1Range" ]
        if ( x2Range is None ): x2Range = config["Data_avgRange"]
    elif ( pltAxis == 2 ):
        pAxiskey = "yAxis"
        meanAxis = 1
        if ( x1Range is None ): x1Range = config["Data_avgRange"]
        if ( x2Range is None ): x2Range = config["Data_x2Range" ]
    Data = fpc.FetchPIC( job =job , jobDir=jobDir, kstep  =kstep, \
                         keys=keys, config=config, x1Range=x1Range, x2Range=x2Range )
    
    ret = { "Axis":Data[pAxiskey] }
    for key in keys:
        ret[key] = np.ravel( np.mean( Data[key], axis=meanAxis ) )
    return( ret )


# ------------------------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    Data    = Fetch1D( job  =args["job"]  , jobDir=args["jobDir"], \
                       kstep=args["kstep"], keys  = ["Bx"] )
    import myStyle.plot1D as pl1
    pl1.plot1D( xAxis=Data["Axis"], yAxis=Data["Bx"] ) 
