import numpy                   as np
import myStyle.LoadConfig      as lcf
import myStyle.plot1D          as pl1
import myStyle.configSettings  as cfs
import myConvert.arr2dct       as a2d
import myBasicAlgs.integrate1D as itg

# ======================================== #
# ===  エリアごとに E.J を評価         === #
# ======================================== #
def AreaEdJAnalysis( job=None, jobDir=None, datDir=None, pngDir=None, pngFile=None, config=None, Flag__integral=True ):
    # ------------------------- #
    # --- [1] 準備          --- #
    # ------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( datDir  is None ): datDir  = "{0}dat/".format( jobDir )
    if ( pngDir  is None ): pngDir  = "{0}png/".format( jobDir )
    if ( pngFile is None ): pngFile = "{0}Area_EdJAnalysis.png".format( pngDir )
    aDict     = {}
    akeys     = ["out0", "out1", "out2", "all", "LCFS"]
    aLegs     = ["$DS$", "$DS_{out}$", "$DS_{in}$","$All$", "LCFS"]
    aLegs     = a2d.arr2dct( Data=aLegs, keys=akeys )
    datFile   = "{0}EdotJAnalysis_{1}.dat"
    with open( datFile.format( datDir, "all" ), "r" ) as f:
        items = ( ( f.readline() ).replace("#","") ).split()
    for akey in akeys:
        Data  = np.transpose( np.loadtxt( datFile.format( datDir, akey ), comments="#" ) )
        Data  = a2d.arr2dct( Data=Data, keys=items )
        if ( Flag__integral ):
            for ikey in items[2:]:
                Data[ikey] = itg.integrate1D( xAxis=Data["ptime"], yAxis=Data[ikey], initial=0.0 )
        aDict[akey] = Data
    tAxis = Data["ptime"]
    # ------------------------- #
    # --- [2] コンフィグ    --- #
    # ------------------------- #
    cfs.configSettings( config=config, configType="plot1D_lateral" )
    cfs.configSettings( config=config, configType="FilterClear"    )
    config["plt_yAutoRange"] = True
    config["plt_yRange"]     = [-1000,4000]
    # ------------------------- #
    # --- [3] プロット      --- #
    # ------------------------- #
    #  --   E.Je            --  #
    FigName = pngFile.replace( ".png", "_e.png" )
    fig     = pl1.plot1D( FigName=FigName, config=config )
    fig.addPlot( xAxis=tAxis, yAxis=( aDict["out0"] )["sJesumE"], label=aLegs["out0"] )
    fig.addPlot( xAxis=tAxis, yAxis=( aDict["out1"] )["sJesumE"], label=aLegs["out1"] )
    fig.addPlot( xAxis=tAxis, yAxis=( aDict["out2"] )["sJesumE"], label=aLegs["out2"] )
    fig.addPlot( xAxis=tAxis, yAxis=( aDict["all" ] )["sJesumE"], label=aLegs["all" ] )
    fig.addPlot( xAxis=tAxis, yAxis=( aDict["LCFS"] )["sJesumE"], label=aLegs["LCFS"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()
    #  --   E.Ji            --  #
    FigName = pngFile.replace( ".png", "_i.png" )
    fig     = pl1.plot1D( FigName=FigName, config=config )
    fig.addPlot( xAxis=tAxis, yAxis=( aDict["out0"] )["sJisumE"], label=aLegs["out0"] )
    fig.addPlot( xAxis=tAxis, yAxis=( aDict["out1"] )["sJisumE"], label=aLegs["out1"] )
    fig.addPlot( xAxis=tAxis, yAxis=( aDict["out2"] )["sJisumE"], label=aLegs["out2"] )
    fig.addPlot( xAxis=tAxis, yAxis=( aDict["all" ] )["sJisumE"], label=aLegs["all" ] )
    fig.addPlot( xAxis=tAxis, yAxis=( aDict["LCFS"] )["sJisumE"], label=aLegs["LCFS"] )
    fig.addLegend()
    fig.setAxis()
    fig.writeFigure()


# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    pngDir  = "./png/AreaEdJAnalysis/"
    ret     = AreaEdJAnalysis( job=args["job"], jobDir=args["jobDir"], pngDir=pngDir )
