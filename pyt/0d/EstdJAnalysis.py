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
import myAnalysis.staticPotential    as spt

# ======================================== #
# ===    J.Eの時間領域 解析            === #
# ======================================== #
def EstdJAnalysis( job    =None, jobDir=None, ksteps =None , CnsFile=None, \
                   datFile=None, datDir=None, pngFile=None , pngDir =None, config=None, \
                   NewFile=True, Flag__CallCalculater=True , Flag__CallPlotter=True, MaskType=None, Flag__integral=True ):
    # ------------------------------- #
    # --- [1]   引数チェック      --- #
    # ------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( ksteps  is None ): ksteps  = np.arange( config["arg_Iter"] )*config["arg_Step"]+config["arg_Init"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( datDir  is None ): datDir  = "{0}dat/".format( jobDir )
    if ( pngDir  is None ): pngDir  = "{0}png/".format( jobDir )
    if ( datFile is None ): datFile = datDir + "EdotJAnalysis_{0}.dat"
    if ( pngFile is None ): pngFile = pngDir + "EdotJAnalysis_{0}.png"
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    datFile_st, datFile_id = datFile.format("st"), datFile.format("id")
    pngFile_st, pngFile_id = pngFile.format("st"), pngFile.format("id")
    # ------------------------------- #
    # --- [2]  J.E を解析         --- #
    # ------------------------------- #
    if ( Flag__CallCalculater ):
        Nsteps = len( ksteps )
        if ( NewFile ):                                # -- 新規ファイル -- #
            items = "# fstep ptime "\
                    "sJepolE sJetorE sJesumE sJeprlE sJepp1E sJepp2E "\
                    "sJipolE sJitorE sJisumE sJiprlE sJipp1E sJipp2E\n"
            with open( datFile_st, "w" ) as f: f.write( items )
            with open( datFile_id, "w" ) as f: f.write( items )
        pgb.progressBar( iteration=0, total=Nsteps )   # -- 進捗表示     -- #
        for i,kstep in enumerate( ksteps ):
            ret = EstdJ_calculator( job=job, jobDir=jobDir, kstep=kstep, config=config, MaskType=MaskType )
            rst = np.array( ret[0:2] + ret[ 2:14] ).reshape( (1,-1) )
            rid = np.array( ret[0:2] + ret[14:26] ).reshape( (1,-1) )
            pgb.progressBar( iteration=i+1, total=Nsteps )
            with open( datFile_st, "ab" ) as f:
                np.savetxt( f, rst )
            with open( datFile_id, "ab" ) as f:
                np.savetxt( f, rid )
    # ------------------------------- #
    # --- [3] 返却 /  描画        --- #
    # ------------------------------- #
    #  -- [3-1] .png 書き出し     --  #
    if ( Flag__CallPlotter ):
        EstdJAnalysis_plotter( datFile=datFile_st, config=config, pngFile=pngFile_st, CnsFile=CnsFile, \
                               Flag__integral=Flag__integral )
        EstdJAnalysis_plotter( datFile=datFile_id, config=config, pngFile=pngFile_id, CnsFile=CnsFile, \
                               Flag__integral=Flag__integral )


# ======================================== #
# ===    J.Eを計算して表示             === #
# ======================================== #
def EstdJ_calculator( job=None, jobDir=None, kstep =None, config=None, MaskType=None, RawData=False ):
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
    keys  = [ "Bx","By","Bz","Jex","Jey","Jez","Jix","Jiy","Jiz","psi" ]
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
    volxMask = volm * mask
    JeppD    = ppD.RefppDecompose( r1=Data["Bx"] , r2=Data["By"] , r3=Data["Bz"] , \
                                   v1=Data["Jex"], v2=Data["Jey"], v3=Data["Jez"]  )
    JippD    = ppD.RefppDecompose( r1=Data["Bx"] , r2=Data["By"] , r3=Data["Bz"] , \
                                   v1=Data["Jix"], v2=Data["Jiy"], v3=Data["Jiz"]  )
    Evec     = spt.staticPotential( job=job, jobDir=jobDir, kstep=kstep, config=config, \
                                    alpha=0.5, nFilter=10, RawData=RawData )
    # ---------------------------- #
    # --- [3] Est.Je ( elec. ) --- #
    # ---------------------------- #
    #  -- [3-1] Je// .E        --  #
    sJeprlEst = np.sum( ( JeppD["vparax"]*Evec["Estx"] + JeppD["vparay"]*Evec["Esty"] + JeppD["vparaz"]*Evec["Estz"] ) * volxMask )
    #  -- [3-2] Je|1 .E        --  #
    sJepp1Est = np.sum( ( JeppD["vprp1x"]*Evec["Estx"] + JeppD["vprp1y"]*Evec["Esty"] + JeppD["vprp1z"]*Evec["Estz"] ) * volxMask )
    #  -- [3-3] Je|2 .E        --  #
    sJepp2Est = np.sum( ( JeppD["vprp2x"]*Evec["Estx"] + JeppD["vprp2y"]*Evec["Esty"] + JeppD["vprp2z"]*Evec["Estz"] ) * volxMask )
    #  -- [3-4] pol,tor,sum    --  #
    sJepolEst= np.sum( ( Data["Jex"]*Evec["Estx"] + Data["Jez"]*Evec["Estz"] ) * volxMask )
    sJetorEst= np.sum( ( Data["Jey"]*Evec["Esty"]                            ) * volxMask )
    sJesumEst= sJepolEst + sJetorEst
    # ---------------------------- #
    # --- [4] Est.Ji (  ion  ) --- #
    # ---------------------------- #
    #  -- [4-1] Ji// .E        --  #
    sJiprlEst = np.sum( ( JippD["vparax"]*Evec["Estx"] + JippD["vparay"]*Evec["Esty"] + JippD["vparaz"]*Evec["Estz"] ) * volxMask )
    #  -- [4-2] Ji|1 .E        --  #
    sJipp1Est = np.sum( ( JippD["vprp1x"]*Evec["Estx"] + JippD["vprp1y"]*Evec["Esty"] + JippD["vprp1z"]*Evec["Estz"] ) * volxMask )
    #  -- [4-3] Ji|2 .E        --  #
    sJipp2Est = np.sum( ( JippD["vprp2x"]*Evec["Estx"] + JippD["vprp2y"]*Evec["Esty"] + JippD["vprp2z"]*Evec["Estz"] ) * volxMask )
    #  -- [4-4] pol,tor,sum    --  #
    sJipolEst= np.sum( ( Data["Jix"]*Evec["Estx"] + Data["Jiz"]*Evec["Estz"] ) * volxMask )
    sJitorEst= np.sum( ( Data["Jiy"]*Evec["Esty"]                            ) * volxMask )
    sJisumEst= sJipolEst + sJitorEst
    # ---------------------------- #
    # --- [3] Eid.Je ( elec. ) --- #
    # ---------------------------- #
    #  -- [3-1] Je// .E        --  #
    sJeprlEid = np.sum( ( JeppD["vparax"]*Evec["Eidx"] + JeppD["vparay"]*Evec["Eidy"] + JeppD["vparaz"]*Evec["Eidz"] ) * volxMask )
    #  -- [3-2] Je|1 .E        --  #
    sJepp1Eid = np.sum( ( JeppD["vprp1x"]*Evec["Eidx"] + JeppD["vprp1y"]*Evec["Eidy"] + JeppD["vprp1z"]*Evec["Eidz"] ) * volxMask )
    #  -- [3-3] Je|2 .E        --  #
    sJepp2Eid = np.sum( ( JeppD["vprp2x"]*Evec["Eidx"] + JeppD["vprp2y"]*Evec["Eidy"] + JeppD["vprp2z"]*Evec["Eidz"] ) * volxMask )
    #  -- [3-4] pol,tor,sum    --  #
    sJepolEid= np.sum( ( Data["Jex"]*Evec["Eidx"] + Data["Jez"]*Evec["Eidz"] ) * volxMask )
    sJetorEid= np.sum( ( Data["Jey"]*Evec["Eidy"]                            ) * volxMask )
    sJesumEid= sJepolEid + sJetorEid
    # ---------------------------- #
    # --- [4] Eid.Ji (  ion  ) --- #
    # ---------------------------- #
    #  -- [4-1] Ji// .E        --  #
    sJiprlEid = np.sum( ( JippD["vparax"]*Evec["Eidx"] + JippD["vparay"]*Evec["Eidy"] + JippD["vparaz"]*Evec["Eidz"] ) * volxMask )
    #  -- [4-2] Ji|1 .E        --  #
    sJipp1Eid = np.sum( ( JippD["vprp1x"]*Evec["Eidx"] + JippD["vprp1y"]*Evec["Eidy"] + JippD["vprp1z"]*Evec["Eidz"] ) * volxMask )
    #  -- [4-3] Ji|2 .E        --  #
    sJipp2Eid = np.sum( ( JippD["vprp2x"]*Evec["Eidx"] + JippD["vprp2y"]*Evec["Eidy"] + JippD["vprp2z"]*Evec["Eidz"] ) * volxMask )
    #  -- [4-4] pol,tor,sum    --  #
    sJipolEid= np.sum( ( Data["Jix"]*Evec["Eidx"] + Data["Jiz"]*Evec["Eidz"] ) * volxMask )
    sJitorEid= np.sum( ( Data["Jiy"]*Evec["Eidy"]                            ) * volxMask )
    sJisumEid= sJipolEid + sJitorEid
    # ---------------------------- #
    # --- [5]     返却         --- #
    # ---------------------------- #
    return( [ fstep    , ptime    , \
              sJepolEst, sJetorEst, sJesumEst, sJeprlEst, sJepp1Est, sJepp2Est, \
              sJipolEst, sJitorEst, sJisumEst, sJiprlEst, sJipp1Est, sJipp2Est, \
              sJepolEid, sJetorEid, sJesumEid, sJeprlEid, sJepp1Eid, sJepp2Eid, \
              sJipolEid, sJitorEid, sJisumEid, sJiprlEid, sJipp1Eid, sJipp2Eid ] )


def EstdJAnalysis_plotter( datFile=None, pngFile=None, config=None, CnsFile=None, Flag__integral=False ):
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
    FigName = pngFile.replace( "EdotJ", "EstdJe_poltor" )
    fig     = pl1.plot1D( FigName=FigName, config=config )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJepolE"], label=Legs["sJepolE"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJetorE"], label=Legs["sJetorE"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJesumE"], label=Legs["sJesumE"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    #  --   E.Ji(pol,tor)   --  #
    config["plt_yRange"]     = [0.0,2000]
    FigName = pngFile.replace( "EdotJ", "EstdJi_poltor" )
    fig     = pl1.plot1D( FigName=FigName, config=config )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJipolE"], label=Legs["sJipolE"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJitorE"], label=Legs["sJitorE"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJisumE"], label=Legs["sJisumE"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    #  --   E.Je(//,|1,|2)   --  #
    config["plt_yRange"]     = [-1000,4000]
    sJetot  = Data["sJeprlE"] + Data["sJepp1E"] + Data["sJepp2E"]
    FigName = pngFile.replace( "EdotJ", "EstdJe_pp" )
    fig     = pl1.plot1D( FigName=FigName, config=config )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJeprlE"], label=Legs["sJeprlE"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJepp1E"], label=Legs["sJepp1E"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJepp2E"], label=Legs["sJepp2E"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=sJetot, label="Total" )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    #  --   E.Ji(//,|1,|2)   --  #
    config["plt_yRange"]     = [0.0,2000]
    sJitot  = Data["sJiprlE"] + Data["sJipp1E"] + Data["sJipp2E"]
    FigName = pngFile.replace( "EdotJ", "EstdJi_pp" )
    fig     = pl1.plot1D( FigName=FigName, config=config )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJiprlE"], label=Legs["sJiprlE"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJipp1E"], label=Legs["sJipp1E"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=Data["sJipp2E"], label=Legs["sJipp2E"] )
    fig.addPlot( xAxis=Data["ptime"], yAxis=sJitot, label="Total" )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()


# ======================================== #
# ===           テスト実行             === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args     = rar.myRecvArgs()
    OutDir   = "./png/EstdJAnalysis/"
    EstdJAnalysis( job   =args["job"]   , jobDir=args["jobDir"], \
                   ksteps=args["ksteps"]  )
