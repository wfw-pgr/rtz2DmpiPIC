import numpy                   as np
import myStyle.LoadConfig      as lcf
import myStyle.StackBarGroup   as sbg
import myStyle.configSettings  as cfs
import myUtils.LoadDatFile     as ldf
import myConvert.pilearr       as pil
import myConvert.arr2dct       as a2d
import mPICsee.TimeUnitConvert as tuc
import myBasicAlgs.integrate1D as itg

# ========================================== #
# ===  エネルギーゲイン比率を棒グラフ化  === #
# ========================================== #
def sEdJdtPlotter( jobBase=None, pngDir=None, pngFile=None, config=None ):
    # ---------------------------------------- #
    # --- [1]   引数チェック               --- #
    # ---------------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( jobBase is None ): jobBase = "../" + config["pic_jobDir"]
    if ( pngDir  is None ): pngDir  = "./png/"
    if ( pngFile is None ): pngFile = "{0}sEdJdtPlot.png".format( pngDir )
    # ------------------ #
    # -- 対象ジョブ名 -- #
    # ------------------ #
    # keys     = ["sJetorE","sJepolE","sJeprlE","sJepp1E","sJepp2E","sJitorE","sJipolE","sJiprlE","sJipp1E","sJipp2E"]
    keys     = ["sJetorE","sJepolE","sJeprlE","sJepp1E","sJitorE","sJipolE","sJiprlE","sJipp1E"]
    jobs     = ["CtrO23_01","CtrO21_01","CtrO22_01","Co-23_01","Co-21_01","Co-22_01"]
    mime     = [ 25, 100, 400, 25, 100, 400 ]
    # datFiles = ["inLCFS_CoHlx.dat", "inLCFS_CtrHlx.dat"]
    # types    = ["Co-"             , "Counter-"         ]

    # ---------------------------------------- #
    # --- [2] データ読込                   --- #
    # ---------------------------------------- #
    #  -- [2-1] データ読込       -- #
    Data  = {}
    for job in jobs:
        jobDir    = "{0}{1}/".format( jobBase, job )
        datFile   = "{0}dat/EdotJAnalysis.dat".format( jobDir )
        hData     = np.transpose( np.loadtxt( datFile, comments="#" ) )
        with open( datFile, "r" ) as f:
            items = ( ( f.readline() ).replace("#","") ).split()
        hData     = a2d.arr2dct( Data=hData, keys=items )
        hData["job"]   = job
        hData["ptime"] = tuc.TimeUnitConvert( kstep=hData["fstep"], unit="wpi", job=job, jobDir=jobDir )
        hData["ctime"] = tuc.TimeUnitConvert( kstep=hData["fstep"], unit="wce", job=job, jobDir=jobDir )
        Data[job] = hData
    #  -- [2-2] データ時系列積分 -- #
    sData = {}
    for job in jobs:
        hData = {}
        rData = Data[job]
        for key in keys:
            hData[key] = itg.integrate1D( xAxis=rData["ctime"], yAxis=rData[key], initial=0.0 )
        sData[job] = hData
    eData = []
    for job in jobs:
        hData = []
        for key in keys:
            hData.append( ( ( sData[job] )[key] )[-1] )
        eData.append( hData )
    # eData = ( ( np.array( eData ) ).reshape( (2,3,-1) ) ).transpose( (1,0,2) )
    eData = ( ( np.array( eData ) ).reshape( (2,3,-1) ) )
    sJetorpolE = eData[:,:,0:2]
    sJeprlppdE = eData[:,:,2:4]
    sJitorpolE = eData[:,:,4:6]
    sJiprlppdE = eData[:,:,6:8]
    tpColor    = ["palevioletred","dodgerblue"]
    # ppColor    = ["gold","darkorchid","aquamarin"]
    ppColor    = ["darkorchid","aquamarine"]
    
    # ---------------------------------------- #
    # --- [3] コンフィグ設定               --- #
    # ---------------------------------------- #
    cfs.configSettings( config=config, configType="FilterClear" )
    Normalize = False
    Noyticks  = False
    
    # ---------------------------------------- #
    # --- [4] プロット                     --- #
    # ---------------------------------------- #
    sJetorpolE = np.abs( sJetorpolE )
    FigName    = pngFile.replace( "Plot.png", "_poltor_e.png" )
    sbg.StackBarGroup( Data=sJetorpolE, FigName=FigName, Normalize=Normalize, colors=tpColor, Noyticks=Noyticks )
    sJeprlppdE = np.abs( sJeprlppdE )
    FigName    = pngFile.replace( "Plot.png", "_prlppd_e.png" )
    sbg.StackBarGroup( Data=sJeprlppdE, FigName=FigName, Normalize=Normalize, colors=ppColor, Noyticks=Noyticks )
    sJitorpolE = np.abs( sJitorpolE )
    FigName    = pngFile.replace( "Plot.png", "_poltor_i.png" )
    sbg.StackBarGroup( Data=sJitorpolE, FigName=FigName, Normalize=Normalize, colors=tpColor, Noyticks=Noyticks )
    sJiprlppdE = np.abs( sJiprlppdE )
    FigName    = pngFile.replace( "Plot.png", "_prlppd_i.png" )
    sbg.StackBarGroup( Data=sJiprlppdE, FigName=FigName, Normalize=Normalize, colors=ppColor, Noyticks=Noyticks )

    # ---------------------------------------- #
    # --- [5] データ百分率プリント         --- #
    # ---------------------------------------- #
    for i in range(2):
        if ( i==0 ): print( "ctr-" )
        if ( i==1 ): print( "co-" )
        for j in range(3):
            if ( j==0 ): print( "mi/me=25" )
            if ( j==1 ): print( "mi/me=100" )
            if ( j==2 ): print( "mi/me=400" )
            print( "[electron]" )
            tpR = sJetorpolE[i,j,:] / np.sum( sJetorpolE[i,j,:] )
            print( "\t toroidal : poloidal      \t {0:>.8f}  : {1:>.8f}  ".format( tpR[0], tpR[1] ) )
            ppR = sJeprlppdE[i,j,:] / np.sum( sJeprlppdE[i,j,:] )
            print( "\t parallel : perpendicular \t {0:>.8f}  : {1:>.8f}  ".format( ppR[0], ppR[1] ) )
            print( "[  ion  ]"  )
            tpR = sJitorpolE[i,j,:] / np.sum( sJitorpolE[i,j,:] )
            print( "\t toroidal : poloidal      \t {0:>.8f}  : {1:>.8f}  ".format( tpR[0], tpR[1] ) )
            ppR = sJiprlppdE[i,j,:] / np.sum( sJiprlppdE[i,j,:] )
            print( "\t parallel : perpendicular \t {0:>.8f}  : {1:>.8f}  ".format( ppR[0], ppR[1] ) )
    
# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    pngDir  = None
    ret     = sEdJdtPlotter( pngDir=pngDir )
