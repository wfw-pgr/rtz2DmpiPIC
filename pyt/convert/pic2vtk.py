import numpy                as np
import myStyle.LoadConfig   as lcf
import mPICsee.retFieldKeys as rfk
import mPICsee.FetchPIC     as fpc
from pyevtk.hl              import gridToVTK

# ======================================== #
# ===    pic => .vtu 用の ルーチン     === #
# ======================================== #
def pic2vtk( job=None, jobDir=None, kstep=None, vtkFile=None, vtkDir=None, keyset=None, config=None, silent=False ):
    # ------------------------- #
    # --- [1] 引数チェック  --- #
    # ------------------------- #
    if ( config  is None ): config   = lcf.LoadConfig()
    if ( job     is None ): job      = config["pic_job"]
    if ( kstep   is None ): kstep    = config["pic_kstep"]
    if ( jobDir  is None ): jobDir   = "{0}{1}/"          .format( config["pic_jobDir"], job )
    if ( vtkDir  is None ): vtkDir   = "{0}vtk/"          .format( jobDir )
    if ( vtkFile is None ): vtkFile  = "{0}picvtk{1:08}".format( vtkDir, kstep  )
    if ( keyset  is None ): keyset   = "Field"
    # ------------------------- #
    # --- [2] データ呼出    --- #
    # ------------------------- #
    #  --    keyset 設定    --  #
    if ( keyset == "Field"     ): Fkeys   = rfk.retFieldKeys( keytype="Field+psi" )
    if ( keyset == "Moment"    ): Fkeys   = rfk.retFieldKeys( keytype="Moment"    )
    if ( keyset == "TAnalysis" ):
        Fkeys   = rfk.retFieldKeys( keytype="Moment"    )
        config["pic_Temperature"] = True
    #  --  ファイル名 設定  --  #
    vtkFile   = vtkFile.replace( "picvtk", keyset )
    #  --    データ読込     --  #
    Data      = fpc.FetchPIC( job=job, jobDir=jobDir, kstep=kstep, keys=Fkeys, config=config )
    #  --    グリッド化     --  #
    xg,yg,zg  = np.meshgrid( Data["xAxis"], Data["yAxis"], np.zeros( (1) ), indexing='ij' )
    pointData = { key :( np.transpose( Data[key] ) )[:,:,np.newaxis] for key in Fkeys }
    # -- ポイントデータとして書き込み -- #
    gridToVTK(vtkFile, xg, yg, zg, pointData = pointData )
    if ( not( silent ) ): print( "[pic2vtk] job= {0}, kstep= {1:8} is processed...".format( args["job"], kstep ) )


# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=='__main__' ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    vtkDir  = None
    if   ( args["key"]    == None        ):
        keyset = None
    elif ( args["key"][0] == "Field"     ):
        keyset = "Field"
    elif ( args["key"][0] == "Moment"    ):
        keyset = "Moment"
    elif ( args["key"][0] == "TAnalysis" ):
        keyset = "TAnaly"
    for kstep in args["ksteps"]:
        vtk = pic2vtk( job=args["job"], jobDir=args["jobDir"], kstep=kstep, vtkDir=vtkDir, keyset=keyset )
