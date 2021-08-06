import numpy                  as np
import myStyle.LoadConfig     as lcf
import myStyle.plot1D         as pl1
import myStyle.configSettings as cfs
import matplotlib.cm          as cm

def potentialWell_plotter( job=None, pngFile=None, datDir=None, config=None, ksteps=None ):
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( datDir  is None ): datDir  = "./potentialWell/"
    if ( pngFile is None ): pngFile = "{0}{1}pWell.png".format( datDir, job )

    cfs.configSettings( config=config, configType="plot1D_def" )
    config["plt_xAutoRange"] = True
    config["plt_yAutoRange"] = False
    config["plt_yRange"]     = [-1.2,0.3]
    fig = pl1.plot1D( config=config, FigName=pngFile )
    for i,kstep in enumerate( ksteps ):
        datFile  = "{0}{1}pWell{2:08}.dat".format( datDir, job, kstep )
        with open( datFile, "r" ) as f:
            Data = np.loadtxt( f )
        fig.addPlot( xAxis=Data[:,0], yAxis=Data[:,1], color=cm.jet(i/90.0) )
    fig.setAxis()
    fig.writeFigure()
    
            
# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    potentialWell_plotter( job=args["job"], ksteps=args["ksteps"] )
