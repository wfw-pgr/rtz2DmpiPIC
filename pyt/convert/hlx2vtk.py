import numpy                as np
import myStyle.LoadConfig   as lcf
import mPICsee.FetchHelix   as fhx
import mPICsee.retFieldKeys as rfk
from pyevtk.hl              import gridToVTK

# -- ヘリシティ を .vts 化するための ルーチン -- #
def hlx2vtk( job=None, jobDir=None, kstep=None, OutFile=None, OutDir=None, config=None, Flag__FieldReturn=True ):
    # ------------------------- #
    # --- [1] 引数チェック  --- #
    # ------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"      .format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}vtk/"      .format( jobDir )
    if ( OutFile is None ): OutFile = "Helix2D{0:08}".format( kstep  )
    OutFile = OutDir + OutFile

    # ------------------------- #
    # --- [2] データ呼出    --- #
    # ------------------------- #
    Data  = fhx.FetchHelix( job   =job   , jobDir=jobDir, kstep =kstep, \
                            config=config, Flag__FieldReturn=Flag__FieldReturn )
    Hkeys = rfk.retFieldKeys( keytype="Helicity+Field" )
    xaxs  = Data["xAxis"]
    yaxs  = Data["yAxis"]
    zaxs  = np.zeros( (1) )
    # -- グリッド化 -- #
    xg,yg,zg  = np.meshgrid( xaxs, yaxs, zaxs, indexing='ij' )
    print( np.transpose( Data["Bx"] )[:,:,np.newaxis].shape )
    pointData = { key :( np.transpose( Data[key] ) )[:,:,np.newaxis] for key in Hkeys }
    # -- ポイントデータとして書き込み -- #
    gridToVTK( OutFile, xg, yg, zg, pointData = pointData )
    
    
#  -----------------------------------------------------------  #
if ( __name__=='__main__' ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    job     = args["job"]
    jobDir  = args["jobDir"]
    OutDir  = "./vts/"
    for kstep in args["ksteps"]:
        print( "[hlx2vtk] job= {0}, kstep= {1:8} is under processing....".format( job, kstep ), end="" )
        vtk = hlx2vtk( job=job, jobDir=jobDir, kstep=kstep, OutDir=OutDir  )
        print( "\t [Done]" )
