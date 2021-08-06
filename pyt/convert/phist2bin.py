import numpy                        as np
import myStyle.LoadConfig           as lcf
import myAnalysis.staticPotential   as spt
import IORoutines.LoadConst         as lcn
import IORoutines.LoadFortranBinary as lfb
import myConvert.arr2dct            as a2d
import myStyle.cMap2D               as clm
import myStyle.configSettings       as cfs

# ================================================================ #
# ===  phist2bin :: 静電ポテンシャルを .bin ファイルへ書込     === #
# ================================================================ #
def phist2bin( job   =None, ksteps =None, jobDir=None, binDir =None, FileName=None, \
               config=None, RawData=True, alpha =None, nFilter=None ):
    # ---------------------------------------- #
    # --- [1]    引数チェック              --- #
    # ---------------------------------------- #
    if ( config   is None ): config   = lcf.LoadConfig()
    if ( job      is None ): job      = config["pic_job"]
    if ( ksteps   is None ): ksteps   = np.arange( config["arg_Iter"] )*config["arg_Step"] + config["arg_Init"]
    if ( jobDir   is None ): jobDir   = "{0}{1}".format( config["pic_jobDir"], job )
    if ( binDir   is None ): binDir   = "{0}bin/".format( jobDir )
    if ( FileName is None ): FileName = "{0}phist__.bin".format( binDir )
    if ( alpha    is None ): alpha    = 0.0
    if ( nFilter  is None ): nFilter  = 0.0
    # ---------------------------------------- #
    # --- [2] 変換.bin 書出し(pointData)   --- #
    # ---------------------------------------- #
    for kstep in ksteps:
        FileName_ = FileName.replace( "phist__","phist{0:08}".format( kstep ) )
        spt.staticPotential( job    =job    , jobDir  =jobDir   , kstep =kstep , alpha=alpha, nFilter=nFilter, \
                             RawData=RawData, FileName=FileName_, config=config, Flag__pointData=False )


# ================================================================ #
# ===  phist2bin2png :: .bin 読み取って .pngへ変換             === #
# ================================================================ #
def phist2bin2png( job=None, jobDir=None, ksteps=None, binDir=None, pngDir=None, \
                   binFile=None, pngFile=None, cnsFile=None, config=None ):
    # ---------------------------------------- #
    # --- [1]    引数チェック              --- #
    # ---------------------------------------- #
    if ( config   is None ): config   = lcf.LoadConfig()
    if ( job      is None ): job      = config["pic_job"]
    if ( jobDir   is None ): jobDir   = "{0}{1}".format( config["pic_jobDir"], job )
    if ( ksteps   is None ): ksteps   = np.arange( config["arg_Iter"] )*config["arg_Step"] + config["arg_Init"]
    if ( binDir   is None ): binDir   = "{0}bin/".format( jobDir )
    if ( pngDir   is None ): pngDir   = "{0}png/".format( jobDir )
    if ( binFile  is None ): binFile  = "{0}phist__.bin".format( binDir )
    if ( pngFile  is None ): pngFile  = "{0}phist__.png".format( pngDir )
    if ( cnsFile  is None ): cnsFile  = "{0}dat/constants.dat".format( jobDir )
    # ---------------------------------------- #
    # --- [2] 変換.bin 書出し(pointData)   --- #
    # ---------------------------------------- #
    const = lcn.LoadConst( FileName=cnsFile )
    keys  = ["phist","Flux","Estx","Esty","Estz","Eidx","Eidy","Eidz"]
    size  = ( 8,const["LJ"],const["LI"] )
    cfs.configSettings( config=config, configType="cMap2D_def" )
    config["cmp_AutoLevel"] = False
    config["cmp_MaxMin"]    = [-0.2,0.2]
    for kstep in ksteps:
        binFile_ = binFile.replace( "phist__","phist{0:08}".format( kstep ) )
        pngFile_ = pngFile.replace( "phist__","phist{0:08}".format( kstep ) )
        rData    = lfb.LoadFortranBinary( FileName=binFile_, size=size )
        Data     = a2d.arr2dct( keys=keys, Data=rData )
        clm.cMap2D( cMap=Data["phist"], config=config, FigName=pngFile_ )


        
# ================================================================ #
# ===  実行部                                                  === #
# ================================================================ #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    phist2bin    ( job=args["job"], jobDir=args["jobDir"], ksteps=args["ksteps"] )
    # phist2bin2png( job=args["job"], jobDir=args["jobDir"], ksteps=args["ksteps"] )
