import numpy              as np
import myStyle.LoadConfig as lcf
import myStyle.pickColor  as pcl


# ================================================================ #
# ===  vdfRegion :: vdf 解析範囲を表示                         === #
# ================================================================ #
def vdfRegion( vcnFile=None, config=None, LI=None, LJ=None, x1Min=None, x1Max=None, x2Min=None, x2Max=None ):
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( vcnFile is None ): vcnFile = "vdfSetup.dat"
    if ( LI      is None ): LI      = config["pic_LI"]
    if ( LJ      is None ): LJ      = config["pic_LJ"]
    if ( x1Min   is None ): x1Min   = 0.0
    if ( x1Max   is None ): x1Max   = float( LI-1 )
    if ( x2Min   is None ): x2Min   = 0.0
    if ( x2Max   is None ): x2Max   = float( LJ-1 )
    x1Leng = x1Max - x1Min
    x2Leng = x2Max - x2Min
    x1Axis = np.linspace( x1Min, x1Max, LI )
    x2Axis = np.linspace( x2Min, x2Max, LJ )
    with open( vcnFile, "r" ) as f:
        line1 = ( f.readline() ).split()
        line2 = ( f.readline() ).split()
        Ndiv  = int( line1[1] )
        Nset  = int( line2[1] )
    vdfcns    = ( np.loadtxt( vcnFile, skiprows=3 ) )
    vdfBox    =  vdfcns[:,0:4]
    vdfBox[:,0:2] = vdfBox[:,0:2] * x1Leng + x1Min
    vdfBox[:,2:4] = vdfBox[:,2:4] * x2Leng + x2Min
    print( vdfBox )
    colors        = pcl.pickColor( nColor=Nset )

    cMap = np.zeros( (LJ,LI) )
    import myStyle.cMap2D as clm
    fig  = clm.cMap2D( xAxis=x2Axis, yAxis=x1Axis, config=config )
    fig.addcMap( cMap=cMap )
    for i in range( Nset ):
        fig.addPlot( xAxis=[vdfBox[i,2],vdfBox[i,2]], yAxis=[vdfBox[i,0],vdfBox[i,1]], color=colors[i] )
    for i in range( Nset ):
        fig.addPlot( xAxis=[vdfBox[i,2],vdfBox[i,3]], yAxis=[vdfBox[i,1],vdfBox[i,1]], color=colors[i] )
    for i in range( Nset ):
        fig.addPlot( xAxis=[vdfBox[i,3],vdfBox[i,3]], yAxis=[vdfBox[i,1],vdfBox[i,0]], color=colors[i] )
    for i in range( Nset ):
        fig.addPlot( xAxis=[vdfBox[i,3],vdfBox[i,2]], yAxis=[vdfBox[i,0],vdfBox[i,0]], color=colors[i] )
    for i in range( Nset ):
        fig.ax1.text( (vdfBox[i,2]+vdfBox[i,3])*0.5, (vdfBox[i,0]+vdfBox[i,1])*0.5, str(i), fontsize=8 )
    fig.writeFigure()


if ( __name__ )==( "__main__"):
    import myUtils.genArgs as gar
    args = gar.genArgs()
    box  = vdfRegion()
