import numpy                as np
import myStyle.LoadConfig   as lcf
import mPICsee.retFieldKeys as rfk
import mPICsee.FetchPIC     as fpc
from pyevtk.hl              import gridToVTK

# -- PIC データを .vtu 化するための ルーチン -- #
def Tei2vtk( job=None, jobDir=None, kstep=None, OutFile=None, OutDir=None, KeySet=None, config=None ):
    # ------------------------- #
    # --- [1] 引数チェック  --- #
    # ------------------------- #
    if ( config  is None ): config   = lcf.LoadConfig()
    if ( job     is None ): job      = config["pic_job"]
    if ( kstep   is None ): kstep    = config["pic_kstep"]
    if ( jobDir  is None ): jobDir   = "{0}{1}/"      .format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir   = "{0}vtk/"      .format( jobDir )
    if ( OutFile is None ): OutFile  = "Field2D{0:08}".format( kstep  )
    if ( KeySet  is None ): KeySet   = "Field"
    if ( KeySet=="Moment"): OutFile  = OutFile.replace( "Field2D", "Moment" )
    OutFile = OutDir + OutFile

    # ------------------------- #
    # --- [2] データ呼出    --- #
    # ------------------------- #
    if ( KeySet == "Moment" ): Fkeys = rfk.retFieldKeys( keytype="Moment" )
    Data      = fpc.FetchPIC( job=job, jobDir=jobDir, kstep=kstep, keys=Fkeys, config=config )
    xaxs      = Data["xAxis"]
    yaxs      = Data["yAxis"]
    zaxs      = np.zeros( (1) )
    Tkeys     = [ "ne0" , "ni0" , \
                  "Texx", "Teyy", "Tezz", "Texy", "Teyz", "Texz", "Tet1", "Ten1", "Ten2", \
                  "Tixx", "Tiyy", "Tizz", "Tixy", "Tiyz", "Tixz", "Tit1", "Tin1", "Tin2"  ]
    pointData         = { "ne0":Data["ne0"], "ni0":Data["ni0"] }
    pointData["Texx"] = Data["pexx"] / 
    
    # -- グリッド化 -- #
    xg,yg,zg  = np.meshgrid( xaxs, yaxs, zaxs, indexing='ij' )
    pointData = { key :( np.transpose( Data[key] ) )[:,:,np.newaxis] for key in Tkeys }    
    # -- ポイントデータとして書き込み -- #
    gridToVTK(OutFile, xg, yg, zg, pointData = pointData )
    
    
#  -----------------------------------------------------------  #
if ( __name__=='__main__' ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    job     = args["job"]
    jobDir  = args["jobDir"]
    OutDir  = None
    if   ( args["key"][0] == "Field"  ):
        KeySet = "Field"
    elif ( args["key"][0] == "Moment" ):
        KeySet = "Moment"
    else:
        KeySet = None
    for kstep in args["ksteps"]:
        print( "[Tei2vtk] job= {0}, kstep= {1:8} is under processing....".format( job, kstep ), end="" )
        vtk = Tei2vtk( job=job, jobDir=jobDir, kstep=kstep, OutDir=OutDir, KeySet=KeySet )
        print( "\t [Done]" )
    
