import numpy                        as np
import myStyle.LoadConfig           as lcf
import myStyle.cMap2D               as clm
import mPICsee.FetchVDF             as fvd
import myStyle.configSettings       as cfs

# ================================================================ #
# ===  FetchVDF :: VDFをロードして返却                         === #
# ================================================================ #
def vdf2cmp( job   =None, kstep =None, jobDir=None, pngDir=None, pngFile=None, config=None, \
             axtype=None, nstype=None, vdfset=None ):
    # ---------------------------------------- #
    # --- [1]    引数チェック              --- #
    # ---------------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( pngDir  is None ): pngDir  = "{0}png/".format( jobDir )
    if ( pngFile is None ): pngFile = "{0}vdf___set___.png".format( pngDir )
    if ( nstype  is None ): nstype  = "e"
    if ( axtype  is None ): axtype  = "yz"
    if ( vdfset  is None ): vdfset  = 2
    # ---------------------------------------- #
    # --- [2] データ読み込み               --- #
    # ---------------------------------------- #
    vdfData  = fvd.FetchVDF( job=job, jobDir=jobDir, kstep=kstep, config=config )
    ns       = { "e" :0, "i" :1         }[nstype]
    ax       = { "xy":0, "yz":1, "zx":2 }[axtype]
    cMap     = ( vdfData["vdf"])[vdfset,ns,ax,:,:]
    print( np.sum(cMap) )
    pngFile_ = pngFile.replace( "vdf___set___","vdf_{0}_{1}_set{2}".format( nstype, axtype, vdfset ) )
    # ---------------------------------------- #
    # --- [3] プロット                     --- #
    # ---------------------------------------- #
    cfs.configSettings( config=config, configType="cMap2D_thesis" )
    clm.cMap2D( cMap=cMap, config=config, FigName=pngFile_ )
    
# ================================================================ #
# ===  実行部                                                  === #
# ================================================================ #
if ( __name__=="__main__" ):
    import myUtils.genArgs as gar
    args    = gar.genArgs()
    vdf2cmp( job   =args["job"] , jobDir=args["jobDir"] , axtype=args["key"], kstep=args["ksteps"][0], \
             nstype=args["mode"], vdfset=args["integer"], pngDir=args["pngDir"] )
