import sys
import numpy                      as np
import nkUtilities.load__config   as lcf
import nkUtilities.cMapTri        as cmt
import nkUtilities.configSettings as cfs


# ========================================================= #
# ===  display                                          === #
# ========================================================= #
def display( datFile=None ):

    x_, y_, z_ = 0, 1, 2
    
    # ------------------------------------------------- #
    # --- [1] Arguments                             --- #
    # ------------------------------------------------- #
    config  = lcf.load__config()
    pngFile = datFile.replace( "dat", "png" )

    # ------------------------------------------------- #
    # --- [2] Fetch Data                            --- #
    # ------------------------------------------------- #
    import nkUtilities.load__pointFile as lpf
    Data  = lpf.load__pointFile( inpFile=datFile, returnType="point" )
    
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
    cmt.cMapTri( xAxis=Data[:,x_], yAxis=Data[:,y_], cMap=Data[:,z_], \
    		 pngFile=pngFile, config=config )

# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):

    import nkUtilities.parse__arguments as arg
    args = arg.parse__arguments()

    if ( args["file"] is None ):
        print( "[display__pointFile.py] --file [datFile name] " )
    else:
        display( datFile=args["file"] )

