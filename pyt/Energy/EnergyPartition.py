import numpy                   as np
import myStyle.LoadConfig      as lcf
import myConvert.pilearr       as pil
import myStyle.StackBarPlot    as sbp
import myStyle.pickColor       as pcl
import mPICsee.TimeUnitConvert as tuc

# ======================================== #
# === 運動エネルギーのプロット         === #
# ======================================== #
def EnergyPartition( job    =None, jobDir=None, ptcFile=None, fldFile     =None, etime=None, \
                     OutFile=None, OutDir=None, config =None, Flag__datOut=True ):
    # -------------------------- #
    # --- [1] 引数チェック   --- #
    # -------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}png/".format( jobDir )
    if ( OutFile is None ): OutFile = "{0}EnergyPartition.png"   .format( OutDir )
    if ( ptcFile is None ): ptcFile = "{0}dat/energy_kinetic.dat".format( jobDir )
    if ( fldFile is None ): fldFile = "{0}dat/energy_field.dat"  .format( jobDir )
    etime = 80.0
    
    # -------------------------- #
    # --- [2] エネルギー読込 --- #
    # -------------------------- #
    #  -- [2-1] データ取得   --  #
    Kenrg  = ( np.loadtxt( ptcFile ) )
    Fenrg  = ( np.loadtxt( fldFile ) )
    if ( etime is None ):
        etime = Fenrg[-1,1]
        etime  = tuc.TimeUnitConvert( Data=etime, job=job, jobDir=jobDir, config=config, unit="wpi" )
    dtime  = tuc.TimeUnitConvert( Data=Kenrg[:,1], job=job, jobDir=jobDir, config=config, unit="wpi" )
    tIndex = np.argmin( np.abs( dtime - etime ) )
    print( dtime.shape, etime, tIndex )
    print( dtime[tIndex], etime )
    Kenrg  = ( np.loadtxt( ptcFile ) )[[1,tIndex],:]
    Fenrg  = ( np.loadtxt( fldFile ) )[[1,tIndex],:]
    Upile  = pil.pilearr( ( Fenrg[:,3:6], Fenrg[:,7:10], Kenrg[:,3:6], Kenrg[:,7:10] ), \
                          NoNewAxis=True, axis=1 )
    Kenrg  = ( np.loadtxt( ptcFile ) )[[1,-1],:]
    Fenrg  = ( np.loadtxt( fldFile ) )[[1,-1],:]
    #  -- [2-2] 規格化       --  #
    for i in range( Upile.shape[0] ):
        Upile[i,:] = Upile[i,:] / np.sum( Upile[i,:] )
    labels = [ "$E_r$"   , "$E_t$"   , "$E_z$"   , "$B_r$"   , "$B_t$"   , "$B_z$"  , \
               "$v_{er}$", "$v_{et}$", "$v_{ez}$", "$v_{ir}$", "$v_{it}$", "$v_{iz}$" ]
    colors = pcl.pickColor( pallete="EBvevi", nColor=12 )
    if ( Flag__datOut ):
        datFile = "{0}EnergyPartition.dat".format( OutDir )
        print( "[EnergyPartition] {0} is saved...".format( datFile ) )
        np.savetxt( datFile, Upile )
    sbp.StackBarPlot( Data   =Upile, xLabel=["Before","After"], labels=labels, colors=colors, \
                      FigName=OutFile )

    
# ======================================== #
# === 実行部                           === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args      = rar.myRecvArgs()
    # OutDir    = "png/EnergyPartition/"
    OutDir    = None
    Output    = EnergyPartition( job=args["job"], jobDir=args["jobDir"], OutDir=OutDir )
