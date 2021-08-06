import numpy                   as np
import myStyle.LoadConfig      as lcf
import myStyle.StackBarGroup   as sbg
import myStyle.configSettings  as cfs
import myUtils.LoadDatFile     as ldf
import myConvert.pilearr       as pil

# ========================================== #
# ===  エネルギーゲイン比率を棒グラフ化  === #
# ========================================== #
def sEdJdtPlotter( pngDir=None, pngFile=None, config=None ):
    # ---------------------------------------- #
    # --- [1]   引数チェック               --- #
    # ---------------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( pngDir  is None ): pngDir  = "./png/"
    if ( pngFile is None ): pngFile = "{0}sEdJdtPlot.png".format( pngDir )
    # ------------------ #
    # -- 対象ジョブ名 -- #
    # ------------------ #
    datFiles = ["inLCFS_CoHlx.dat", "inLCFS_CtrHlx.dat"]
    # types    = ["Co-"             , "Counter-"         ]

    # ---------------------------------------- #
    # --- [2] データ読込                   --- #
    # ---------------------------------------- #
    Data        = {}
    Data["Co-"] = ldf.LoadDatFile( InpFile=datFiles[0] )
    Data["Ctr"] = ldf.LoadDatFile( InpFile=datFiles[1] )
    print( Data["Co-"]["sJetorE"] )
    print( Data["Ctr"]["sJetorE"] )
    torpol_e_CoH = pil.pilearr( ( ( Data["Co-"] )["sJetorE"], ( Data["Co-"] )["sJepolE"] ) )
    torpol_e_Ctr = pil.pilearr( ( ( Data["Ctr"] )["sJetorE"], ( Data["Ctr"] )["sJepolE"] ) )
    torpol_e     = pil.pilearr( ( torpol_e_CoH, torpol_e_Ctr ) ).transpose( (2,0,1) )
    torpol_i_CoH = pil.pilearr( ( ( Data["Co-"] )["sJitorE"], ( Data["Co-"] )["sJipolE"] ) )
    torpol_i_Ctr = pil.pilearr( ( ( Data["Ctr"] )["sJitorE"], ( Data["Ctr"] )["sJipolE"] ) )
    torpol_i     = pil.pilearr( ( torpol_i_CoH, torpol_i_Ctr ) ).transpose( (2,0,1) )
    
    # ---------------------------------------- #
    # --- [3] コンフィグ設定               --- #
    # ---------------------------------------- #
    cfs.configSettings( config=config, configType="FilterClear" )

    # ---------------------------------------- #
    # --- [4] プロット                     --- #
    # ---------------------------------------- #
    torpol_e = np.abs( torpol_e )
    FigName  = pngFile.replace( "Plot.png", "_e.png" )
    sbg.StackBarGroup( Data=torpol_e, FigName=FigName, Normalize=False )
    torpol_i = np.abs( torpol_i )
    FigName  = pngFile.replace( "Plot.png", "_i.png" )
    sbg.StackBarGroup( Data=torpol_i, FigName=FigName, Normalize=False )
    
# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    pngDir  = None
    ret     = sEdJdtPlotter( pngDir=pngDir )
