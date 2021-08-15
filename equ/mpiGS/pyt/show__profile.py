import sys
import numpy                      as np
import nkUtilities.load__config   as lcf
import nkUtilities.cMapTri        as cmt
import nkUtilities.configSettings as cfs


# ========================================================= #
# ===  display                                          === #
# ========================================================= #
def display():

    x_,y_,z_ = 0, 1, 2
    
    # ------------------------------------------------- #
    # --- [1] Arguments                             --- #
    # ------------------------------------------------- #
    import nkUtilities.genArgs as gar
    args    = gar.genArgs()
    job     = args["job"]
    mode    = args["mode"]
    config  = lcf.load__config()
    if ( job  is None ): sys.exit( "[show__profile.py] please specify job  with --job option. " )
    if ( mode is None ): mode = "BgD"
    
    # ------------------------------------------------- #
    # --- [2] Fetch Data                            --- #
    # ------------------------------------------------- #
    datDir  = "job/{0}/{1}".format( job, mode )
    bfdFile = "{0}/Bfd.bin".format( datDir )
    jcrFile = "{0}/Jcr.bin".format( datDir )
    frpFile = "{0}/frp.bin".format( datDir )
    uvcFile = "{0}/uvc.bin".format( datDir )

    import nkUtilities.load__constants as lcn
    prmFile = "{0}/parameter.dat".format( datDir )
    const   = lcn.load__constants( inpFile=prmFile )
    x1a     = np.linspace( const["x1Min"], const["x1Max"], const["LI"] )
    x2a     = np.linspace( const["x2Min"], const["x2Max"], const["LJ"] )
    x1g,x2g = np.meshgrid( x1a, x2a )
    x1g,x2g = np.reshape( x1g, (-1,) ), np.reshape( x2g, (-1,) )
    shape   = (const["LJ"],const["LI"],3)
    
    import nkUtilities.load__fortranBinary as lfb
    bfd     = lfb.load__fortranBinary( inpFile=bfdFile, shape=shape )
    jcr     = lfb.load__fortranBinary( inpFile=jcrFile, shape=shape )
    frp     = lfb.load__fortranBinary( inpFile=frpFile, shape=shape )
    uvc     = lfb.load__fortranBinary( inpFile=uvcFile, shape=shape )
    bfd     = np.reshape( bfd, (-1,3) )
    jcr     = np.reshape( jcr, (-1,3) )
    frp     = np.reshape( frp, (-1,3) )
    uvc     = np.reshape( uvc, (-1,3) )
    
    # ------------------------------------------------- #
    # --- [3] config Settings                       --- #
    # ------------------------------------------------- #
    cfs.configSettings( configType="cMap_def", config=config )
    config["FigSize"]        = (5,5)
    config["cmp_position"]   = [0.16,0.12,0.97,0.88]
    config["xTitle"]         = "X (m)"
    config["yTitle"]         = "Y (m)"
    config["cmp_xAutoRange"] = True
    config["cmp_yAutoRange"] = True
    config["cmp_xRange"]     = [-5.0,+5.0]
    config["cmp_yRange"]     = [-5.0,+5.0]

    # ------------------------------------------------- #
    # --- [4] plot Figure                           --- #
    # ------------------------------------------------- #
    pngFile = datDir + "/bfd_{0}.png"
    cmt.cMapTri( xAxis=x1g, yAxis=x2g, cMap=bfd[:,x_], \
                 pngFile=pngFile.format( "x" ), config=config )
    cmt.cMapTri( xAxis=x1g, yAxis=x2g, cMap=bfd[:,y_], \
                 pngFile=pngFile.format( "y" ), config=config )
    cmt.cMapTri( xAxis=x1g, yAxis=x2g, cMap=bfd[:,z_], \
                 pngFile=pngFile.format( "z" ), config=config )

# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    display()

