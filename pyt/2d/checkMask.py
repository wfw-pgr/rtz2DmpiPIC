import mPICsee.FetchPIC       as fpc
import myStyle.LoadConfig     as lcf
import myStyle.cMap2D         as clm
import myStyle.configSettings as cfs
import myAnalysis.createMask  as msk
import myAnalysis.outletMask  as olt

# ======================================== #
# ===    J.Eを計算して表示             === #
# ======================================== #
def checkMask( job=None, jobDir =None, kstep=None, config=None, pngDir=None ):
    # ------------------------------- #
    # --- [1] 引数チェック        --- #
    # ------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( pngDir  is None ): pngDir  = "{0}png/".format( jobDir )
    # ------------------------------- #
    # --- [2]  データ読込 / 描画  --- #
    # ------------------------------- #
    #  -- [2-1] .png 書き出し     --  #
    Data  = fpc.FetchPIC( job =job    , jobDir=jobDir, kstep =kstep, \
                          keys=["psi"], config=config )
    Flux  = Data["psi"]
    masks = msk.createMask( Flux=Flux, Flag__privateFluxMode=True )
    omask = olt.outletMask( Flux=Flux, psiSurpass=0.1             )
    cfs.configSettings( config=config, configType="cMap2D_def" )
    clm.cMap2D( cMap=masks["fMask0"] , Cntr=Flux, FigName=pngDir+"fMask0__{0:08}.png".format( kstep ) )
    clm.cMap2D( cMap=masks["fMask1"] , Cntr=Flux, FigName=pngDir+"fMask1__{0:08}.png".format( kstep ) )
    clm.cMap2D( cMap=masks["fMask2"] , Cntr=Flux, FigName=pngDir+"fMask2__{0:08}.png".format( kstep ) )
    clm.cMap2D( cMap=masks["fMask3"] , Cntr=Flux, FigName=pngDir+"fMask3__{0:08}.png".format( kstep ) )
    clm.cMap2D( cMap=omask["outlet1"], Cntr=Flux, FigName=pngDir+"outlet1_{0:08}.png".format( kstep ) )
    clm.cMap2D( cMap=omask["outlet2"], Cntr=Flux, FigName=pngDir+"outlet2_{0:08}.png".format( kstep ) )
    return( { "masks":masks, "omask":omask } )

# ======================================== #
# ===           テスト実行             === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    pngDir  = "./checkMask/"
    checkMask( job  =args["job"]  , jobDir=args["jobDir"], \
               kstep=args["kstep"], pngDir=pngDir )
