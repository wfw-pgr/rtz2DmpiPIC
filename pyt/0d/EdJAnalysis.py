import numpy                         as np
import myStyle.plot1D                as pl1
import mPICsee.FetchPIC              as fpc
import myStyle.LoadConfig            as lcf
import myStyle.configSettings        as cfs
import myBasicAlgs.cylindricalVolume as cvl
import myUtils.progressBar           as pgb
import mPICsee.TimeUnitConvert       as tuc
import myConvert.arr2dct             as a2d
import myAnalysis.outletMask         as olm
import myBasicAlgs.integrate1D       as itg
import fLIB.fLIB__RefppDecompose     as ppD
import myAnalysis.createMask         as msk

# ======================================== #
# ===    J.Eの時間領域 解析            === #
# ======================================== #
def EdJAnalysis( job    =None, jobDir=None, ksteps =None , CnsFile=None, \
                 datFile=None, datDir=None, pngFile=None , pngDir =None, config=None, \
                 NewFile=True, Flag__CallCalculater=True , Flag__CallPlotter=True, MaskType="LCFS", Flag__integral=True, RawData=True ):
    # ------------------------------- #
    # --- [1]   引数チェック      --- #
    # ------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( ksteps  is None ): ksteps  = np.arange( config["arg_Iter"] )*config["arg_Step"]+config["arg_Init"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( datDir  is None ): datDir  = "{0}dat/".format( jobDir )
    if ( pngDir  is None ): pngDir  = "{0}png/".format( jobDir )
    if ( datFile is None ): datFile = datDir + "EdotJAnalysis.dat"
    if ( pngFile is None ): pngFile = pngDir + "EdotJAnalysis.png"
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    # ------------------------------- #
    # --- [2]  J.E を解析         --- #
    # ------------------------------- #
    if ( Flag__CallCalculater ):
        Nsteps = len( ksteps )
        if ( NewFile ):                                # -- 新規ファイル -- #
            items = "# fstep ptime "\
                    "sJepolE sJetorE sJesumE sJeprlE sJepp1E sJepp2E "\
                    "sJipolE sJitorE sJisumE sJiprlE sJipp1E sJipp2E\n"
            with open( datFile, "w" ) as f: f.write( items )
        pgb.progressBar( iteration=0, total=Nsteps )   # -- 進捗表示     -- #
        for i,kstep in enumerate( ksteps ):
            ret = EdJ_calculator( job=job, jobDir=jobDir, kstep=kstep, config=config, MaskType=MaskType )
            pgb.progressBar( iteration=i+1, total=Nsteps )
            with open( datFile, "ab" ) as f:
                np.savetxt( f, np.array( ret ).reshape( (1,-1) ) )
    # ------------------------------- #
    # --- [3] 返却 /  描画        --- #
    # ------------------------------- #
    #  -- [3-1] .png 書き出し     --  #
    if ( Flag__CallPlotter ):
        EdJAnalysis_plotter( datFile=datFile, config=config, pngFile=pngFile, CnsFile=CnsFile, \
                             Flag__integral=Flag__integral )


# ======================================== #
# ===    J.Eを計算して表示             === #
# ======================================== #
def EdJ_calculator( job  =None, jobDir=None, kstep =None, config=None, MaskType=None, RawData=False ):
    # -------------------------- #
    # --- [1]  引数チェック  --- #
    # -------------------------- #
    if ( config  is None ): config = lcf.LoadConfig()
    if ( job     is None ): job    = config["pic_job"]
    if ( kstep   is None ): kstep  = config["pic_kstep"]
    if ( jobDir  is None ): jobDir = "{0}{1}/".format( config["pic_jobDir"], job )
    # -------------------------- #
    # --- [2] データ呼び出し --- #
    # -------------------------- #
    #  -- [2-1] 読込         --  #
    fstep = float( kstep )
    ptime = tuc.TimeUnitConvert( job=job, jobDir=jobDir, kstep=kstep, unit="wpi" )
    keys  = [ "Bx","By","Bz","Ex","Ey","Ez","Jex","Jey","Jez","Jix","Jiy","Jiz","psi" ]
    config["cmp_nFilter"]    = 10
    config["cmp_LinearFilt"] = 0.2
    Data  = fpc.FetchPIC( job=job , jobDir=jobDir, kstep =kstep, keys=keys, \
                          FilterKeys=keys, RawData=RawData, config=config )
    xAxis = Data["xAxis"]
    yAxis = Data["yAxis"]
    Flux  = Data["psi"  ]
    #  -- [2-2] 準備         --  #
    if   ( MaskType == "outlet0" ):
        mask = ( olm.outletMask( Flux=Flux ) )["outlet0"]
    elif ( MaskType == "outlet1" ):
        mask = ( olm.outletMask( Flux=Flux ) )["outlet1"]
    elif ( MaskType == "outlet2" ):
        mask = ( olm.outletMask( Flux=Flux ) )["outlet2"]
    elif ( MaskType == "LCFS"    ):
        mask = msk.createMask( Flux=Flux, Flag__LCFSMode=True )
    else:
        mask = 1.0
    volm     = cvl.cylindricalVolume( zAxis=xAxis, rAxis=yAxis )
    JeppD    = ppD.RefppDecompose( r1=Data["Bx"] , r2=Data["By"] , r3=Data["Bz"] , \
                                   v1=Data["Jex"], v2=Data["Jey"], v3=Data["Jez"]  )
    JippD    = ppD.RefppDecompose( r1=Data["Bx"] , r2=Data["By"] , r3=Data["Bz"] , \
                                   v1=Data["Jix"], v2=Data["Jiy"], v3=Data["Jiz"]  )
    # ---------------------------- #
    # --- [3]  E.Je ( elec. )  --- #
    # ---------------------------- #
    #  -- [3-1] Je// .E        --  #
    JeprlE   = JeppD["vparax"]*Data["Ex"] + JeppD["vparay"]*Data["Ey"] + JeppD["vparaz"]*Data["Ez"]
    sJeprlE  = np.sum( JeprlE * volm * mask )
    #  -- [3-2] Je|1 .E        --  #
    Jepp1E   = JeppD["vprp1x"]*Data["Ex"] + JeppD["vprp1y"]*Data["Ey"] + JeppD["vprp1z"]*Data["Ez"]
    sJepp1E  = np.sum( Jepp1E * volm * mask )
    #  -- [3-3] Je|2 .E        --  #
    Jepp2E   = JeppD["vprp2x"]*Data["Ex"] + JeppD["vprp2y"]*Data["Ey"] + JeppD["vprp2z"]*Data["Ez"]
    sJepp2E  = np.sum( Jepp2E * volm * mask )
    #  -- [3-4] pol,tor,sum    --  #
    JepolE   = Data["Jex"]*Data["Ex"] + Data["Jez"]*Data["Ez"]
    JetorE   = Data["Jey"]*Data["Ey"]
    JesumE   = JepolE + JetorE
    sJepolE  = np.sum( JepolE * volm * mask )
    sJetorE  = np.sum( JetorE * volm * mask )
    sJesumE  = np.sum( JesumE * volm * mask )
    # ---------------------------- #
    # --- [4]  E.Ji (  ion  )  --- #
    # ---------------------------- #
    #  -- [4-1] Ji// .E        --  #
    JiprlE   = JippD["vparax"]*Data["Ex"] + JippD["vparay"]*Data["Ey"] + JippD["vparaz"]*Data["Ez"]
    sJiprlE  = np.sum( JiprlE * volm * mask )
    #  -- [4-2] Ji|1 .E        --  #
    Jipp1E   = JippD["vprp1x"]*Data["Ex"] + JippD["vprp1y"]*Data["Ey"] + JippD["vprp1z"]*Data["Ez"]
    sJipp1E  = np.sum( Jipp1E * volm * mask )
    #  -- [4-3] Ji|2 .E        --  #
    Jipp2E   = JippD["vprp2x"]*Data["Ex"] + JippD["vprp2y"]*Data["Ey"] + JippD["vprp2z"]*Data["Ez"]
    sJipp2E  = np.sum( Jipp2E * volm * mask )
    #  -- [4-4] pol,tor,sum    --  #
    JipolE   = Data["Jix"]*Data["Ex"] + Data["Jiz"]*Data["Ez"]
    JitorE   = Data["Jiy"]*Data["Ey"]
    JisumE   = JipolE + JitorE
    sJipolE  = np.sum( JipolE * volm * mask )
    sJitorE  = np.sum( JitorE * volm * mask )
    sJisumE  = np.sum( JisumE * volm * mask )
    # ---------------------------- #
    # --- [5]     返却         --- #
    # ---------------------------- #
    return( [ fstep  , ptime, \
              sJepolE, sJetorE, sJesumE, sJeprlE, sJepp1E, sJepp2E, \
              sJipolE, sJitorE, sJisumE, sJiprlE, sJipp1E, sJipp2E ] )


def EdJAnalysis_plotter( datFile=None, pngFile=None, config=None, CnsFile=None, Flag__integral=False ):
    # ------------------------- #
    # --- [1] 準備          --- #
    # ------------------------- #
    Data = np.transpose( np.loadtxt( datFile, comments="#" ) )
    with open( datFile, "r" ) as f:
        items = ( ( f.readline() ).replace("#","") ).split()
    Legs = [ "fstep","Time", \
             "$(J_e \cdot E)_{pol}$" , "$(J_e \cdot E)_{tor}$"  , "$(J_e \cdot E)_{tot}$"  ,\
             "$(J_e \cdot E)_{para}$", "$(J_e \cdot E)_{perp1}$", "$(J_e \cdot E)_{perp2}$",\
             "$(J_i \cdot E)_{pol}$" , "$(J_i \cdot E)_{tor}$"  , "$(J_i \cdot E)_{tot}$"  ,\
             "$(J_i \cdot E)_{para}$", "$(J_i \cdot E)_{perp1}$", "$(J_i \cdot E)_{perp2}$" ]
    Data = a2d.arr2dct( Data=Data, keys=items )
    Legs = a2d.arr2dct( Data=Legs, keys=items )
    if ( Flag__integral ):
        for key in items[2:]: Data[key] = itg.integrate1D( xAxis=Data["ptime"], yAxis=Data[key], initial=0.0 )
    cfs.configSettings( config=config, configType="plot1D_lateral" )
    cfs.configSettings( config=config, configType="FilterClear"    )
    config["plt_yAutoRange"] = True
    # ------------------------- #
    # --- [2] プロット      --- #
    # ------------------------- #
    #  --   E.Je(pol,tor)   --  #
    config["plt_yRange"]     = [-1000,4000]
    FigName = pngFile.replace( "EdotJ", "EdJe_poltor" )
    fig     = pl1.plot1D( FigName=FigName, config=config )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJepolE"], label=Legs["sJepolE"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJetorE"], label=Legs["sJetorE"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJesumE"], label=Legs["sJesumE"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    print( "[time_EdJ] {0} is outputed...".format( FigName ) )
    #  --   E.Ji(pol,tor)   --  #
    config["plt_yRange"]     = [0.0,2000]
    FigName = pngFile.replace( "EdotJ", "EdJi_poltor" )
    fig     = pl1.plot1D( FigName=FigName, config=config )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJipolE"], label=Legs["sJipolE"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJitorE"], label=Legs["sJitorE"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJisumE"], label=Legs["sJisumE"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    print( "[time_EdJ] {0} is outputed...".format( FigName ) )
    #  --   E.Je(//,|1,|2)   --  #
    config["plt_yRange"]     = [-1000,4000]
    sJetot  = Data["sJeprlE"] + Data["sJepp1E"] + Data["sJepp2E"]
    FigName = pngFile.replace( "EdotJ", "EdJe_pp" )
    fig     = pl1.plot1D( FigName=FigName, config=config )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJeprlE"], label=Legs["sJeprlE"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJepp1E"], label=Legs["sJepp1E"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJepp2E"], label=Legs["sJepp2E"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=sJetot, label="Total" )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    print( "[time_EdJ] {0} is outputed...".format( FigName ) )
    #  --   E.Ji(//,|1,|2)   --  #
    config["plt_yRange"]     = [0.0,2000]
    sJitot  = Data["sJiprlE"] + Data["sJipp1E"] + Data["sJipp2E"]
    FigName = pngFile.replace( "EdotJ", "EdJi_pp" )
    fig     = pl1.plot1D( FigName=FigName, config=config )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJiprlE"], label=Legs["sJiprlE"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJipp1E"], label=Legs["sJipp1E"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJipp2E"], label=Legs["sJipp2E"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=sJitot, label="Total" )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    print( "[time_EdJ] {0} is outputed...".format( FigName ) )


# ======================================== #
# ===           テスト実行             === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args     = rar.myRecvArgs()
    OutDir   = "./png/EdJAnalysis/"
    EdJAnalysis( job   =args["job"]   , jobDir=args["jobDir"], \
                 ksteps=args["ksteps"]  )
