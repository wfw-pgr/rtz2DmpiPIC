import sys
import numpy                  as np
import myStyle.LoadConfig     as lcf
import myStyle.plot1D         as pl1
import myStyle.configSettings as cfs
import myStyle.configRangeSet as crs

def pDiff( config=None ):
    if ( config  is None ): config  = lcf.LoadConfig()
    # ---------------------------------------- #
    # --- [1] データ 読込                  --- #
    # ---------------------------------------- #
    datFile1  = "potentialDiff_CtrO22_01.dat"
    datFile2  = "potentialDiff_CtrO21_01.dat"
    datFile3  = "potentialDiff_CtrO23_01.dat"
    datFiles  = [ datFile1, datFile2, datFile3 ]
    rData     = []
    for datFile in datFiles:
        with open( datFile, "r" ) as f:
            rData.append( ( np.loadtxt( f ) ).transpose() )
    
    # ---------------------------------------- #
    # --- [2] 電位差 / 電場                --- #
    # ---------------------------------------- #
    mRate      = 0.45
    pData      = np.zeros( (3,len(datFiles) ) )
    pData[0,:] = np.array( [ 25.0, 100.0, 400.0 ] )
    for i,rD in enumerate( rData ):
        idx        = np.argmin( np.abs( rD[2,:] - mRate ) )
        pData[1,i] = rD[13,idx]
        pData[2,i] = rD[14,idx]
        print( idx, rD[2,idx], rD[13,idx], rD[14,idx] )
    pData[2,:] = pData[2,:] * np.sqrt( pData[0,:] )
        
    cfs.configSettings( config=config, configType="plot1D_def"  )
    cfs.configSettings( config=config, configType="FilterClear" )
    cfs.configSettings( config=config, configType="plot1D_mark" )
    crs.configRangeSet( config=config, xRange=[0,500], yRange=[0.0,1.0] )
    config["MinimalOut"] = True
    config["xTitle"]     = ""
    config["yTitle"]     = ""
    config["xMajor_Nticks"] = 6
    # config["xMajor_Nticks"] = 6
    # -- -- #
    config["plt_color"]  = "Navy"
    config["plt_marker"] = "o"
    FigName   = "Vdif.png"
    pl1.plot1D( xAxis=pData[0,:], yAxis=pData[1,:], config=config, FigName=FigName )
    # -- -- #
    config["plt_color"]  = "Crimson"
    config["plt_marker"] = "D"
    crs.configRangeSet( config=config, xRange=[0,500], yRange=[0.0,2.0] )
    FigName   = "Eeff.png"
    pl1.plot1D( xAxis=pData[0,:], yAxis=pData[2,:], config=config, FigName=FigName )

# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    OutDir  = "./"
    ret     = pDiff()
