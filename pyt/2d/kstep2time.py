import sys
import numpy              as np
import myUtils.LoadConst  as lcn
import myStyle.LoadConfig as lcf

def TimeUnitConvert( kstep  =None, CnsFile=None, job=None, jobDir=None, config=None, unit=None, display=True ):
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"             .format( config["pic_jobDir"], job )
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    if ( kstep   is None ): sys.exit( "[kstep2time] kstep = ?? " )
    if ( unit    is None ): unit    = config["phys_TimeUnit"]
    Const = lcn.LoadConst( InpFile=CnsFile )
    coef = 1.0
    if ( unit == "wce" ): coef = Const["dt"] 
    if ( unit == "wci" ): coef = Const["dt"] /   Const["mr"]
    if ( unit == "wpe" ): coef = Const["dt"] *   Const["wpewce"]
    if ( unit == "wpi" ): coef = Const["dt"] *   Const["wpewce"] / np.sqrt( Const["mr"] )
    if ( unit == "tAe" ): coef = Const["dt"] / ( Const["wpewce"] * Const["x2Leng"] ) 
    if ( unit == "tAi" ): coef = Const["dt"] / ( Const["wpewce"] * Const["x2Leng"] * np.sqrt( Const["mr"] ) )
    time = kstep * coef
    if ( display ): print ( "[TimeUnitConvert] kstep = {0} :: t = {1} :: ( unit = {2} )".format( kstep, time, unit ) )
    return( time )


if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args = rar.myRecvArgs()
    unit = None
    if ( args["key"][0] == "wce" ): unit = "wce"
    if ( args["key"][0] == "wci" ): unit = "wci"
    if ( args["key"][0] == "wpe" ): unit = "wpe"
    if ( args["key"][0] == "wpi" ): unit = "wpi"
    if ( args["key"][0] == "tAe" ): unit = "tAe"
    if ( args["key"][0] == "tAi" ): unit = "tAi"
    TimeUnitConvert( kstep=args["kstep"], job=args["job"], jobDir=args["jobDir"], unit=unit )
