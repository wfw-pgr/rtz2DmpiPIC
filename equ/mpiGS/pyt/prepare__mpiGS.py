import os, sys
import numpy as np
import nkUtilities.load__constants as lcn
import nkUtilities.save__namelist  as nml

# ========================================================= #
# ===  prepare__mpiGS.py                                === #
# ========================================================= #

def prepare__mpiGS():

    inpFile = "dat/mpiGS.conf"
    nmlFile = "dat/mpiGS.lst"
    
    # ------------------------------------------------- #
    # --- [1] load mpiGS.conf                       --- #
    # ------------------------------------------------- #
    const   = lcn.load__constants( inpFile=inpFile )

    # ------------------------------------------------- #
    # --- [2] save as fortran's namelist file       --- #
    # ------------------------------------------------- #
    groups  = [ "equSettings", "grid", "pic", "equ" ]
    ret     = nml.save__namelist( outFile=nmlFile, const=const, groups=groups )

    # ------------------------------------------------- #
    # --- [3] make job directory                    --- #
    # ------------------------------------------------- #
    jobDir  = "job/{}".format( const["equSettings.job"] )
    os.makedirs( jobDir                        , exist_ok=True )
    os.makedirs( os.path.join( jobDir, "dat/" ), exist_ok=True )
    os.makedirs( os.path.join( jobDir, "png/" ), exist_ok=True )
    os.makedirs( os.path.join( jobDir, "src/" ), exist_ok=True )
    os.makedirs( os.path.join( jobDir, "BgD/" ), exist_ok=True )
    os.makedirs( os.path.join( jobDir, "RgD/" ), exist_ok=True )
    
    return()



# ========================================================= #
# ===   Execution of Pragram                            === #
# ========================================================= #
if ( __name__=="__main__" ):
    prepare__mpiGS()

