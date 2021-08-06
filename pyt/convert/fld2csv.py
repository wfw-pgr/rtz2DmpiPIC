import numpy              as np
import myConvert.bin2csv  as b2c
import myStyle.LoadConfig as lcf

def fld2csv( job   =None, keys  =None, jobDir =None, kstep=None, \
             InpDir=None, OutDir=None, PEiFile=None, csv  =True, \
             nVar  = 18 , config=None ):
    
    # --- [1] 引数設定 --- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "../job/{0}".format( job )
    if ( keys    is None ): keys    = [ "Bx", "By", "Bz" ]
    if ( InpDir  is None ): InpDir  = "{0}/bin/kstep{1:08}/" .format( jobDir, kstep )
    if ( OutDir  is None ): OutDir  = "./csv/"
    if ( PEiFile is None ): PEiFile = "{0}/dat/PEinfo.dat"   .format( jobDir )
    
    # --- [2] キー設定 --- #
    keyDict = { "Bx" :0 , "By" :1 , "Bz" :2 ,
                "Ex" :3 , "Ey" :4 , "Ez" :5 ,
                "Jx" :6 , "Jy" :7 , "Jz" :8 ,
                "Jex":9 , "Jey":10, "Jez":11,
                "Jix":12, "Jiy":13, "Jiz":14,
                "ne" :15, "ni" :16, "qn" :17 }
    index   = [ keyDict[key] for key in keys ]
    PEinfo  = np.array( np.loadtxt( PEiFile ), dtype=np.int32 )
    PEtot   = PEinfo.shape[0]
    # --- [3] bin ==> csv --- #
    for rank in range( PEtot ):
        InpFile  = "Field2D{0:08}_{1:06}.bin".format( kstep, rank )
        size     = [ PEinfo[rank,3], PEinfo[rank,4], nVar ]
        keystack = "("
        for key in keys: keystack = keystack + " {0} ".format(key)
        keystack = keystack + ")"
        print( "[fld2csv] {0} -- ( {1} )  is under processing...".format(InpFile,keystack), end="" )
        b2c.bin2csv( InpFile=InpFile, InpDir=InpDir, OutDir=OutDir , \
                     index  =index  , csv   =csv   , size  =size   , config=config )
        print( "\t [Done]" )

if ( __name__=='__main__' ):
    import myUtils.myRecvArgs as rar
    args = rar.myRecvArgs()
    exe  = fld2csv( job=args["job"], keys=args["key"], kstep=args["ksteps"][0] )
