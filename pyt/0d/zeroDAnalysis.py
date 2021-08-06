import timeCurrent
import timeFlux
import timeMomentum
import myStyle.LoadConfig as lcf

def zeroDAnalysis( job=None, jobDir=None, OutDir=None, config=None ):
    # ----------------- #
    # -- 引数チェック -- #
    # ----------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}png/".format( jobDir )
    # ----------------- #
    # -- ルーチン呼ぶ -- #
    # ----------------- #
    timeCurrent .timeCurrent ( job=job, jobDir=jobDir, OutDir=OutDir )
    print( "[timeCurrent]" )
    timeFlux    .timeFlux    ( job=job, jobDir=jobDir, OutDir=OutDir )
    print( "[timeFlux]" )
    timeMomentum.timeMomentum( job=job, jobDir=jobDir, OutDir=OutDir )
    print( "[timeMomentum]" )

if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    OutDir = "./png/"
    zeroDAnalysis( job=args["job"], OutDir=OutDir )
