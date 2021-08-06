import mPICsee.FetchPIC    as fpc
import myStyle.LoadConfig  as lcf
import myConvert.pilearr   as pil
import fLIB.fLIB__helicity as ghl
import fLIB.fLIB__curl     as crl

def FetchHelix( job  =None, jobDir=None, kstep=None, config=None, Flag__FieldReturn=False ):
    # -------------------------- #
    # --- [1]  引数チェック  --- #
    # -------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/"    .format( config["pic_jobDir"], job )
    
    # -------------------------- #
    # --- [2] データ呼び出し --- #
    # -------------------------- #
    #  -- [2-1]      場      --  #
    Fkeys   = ["Bx","By","Bz","uix","uiy","uiz","uex","uey","uez","psi"]
    Data    = fpc.FetchPIC( job =job  , jobDir=jobDir, kstep =kstep, \
                            keys=Fkeys, config=config )
    zAxis   = Data["xAxis"]
    rAxis   = Data["yAxis"]
    flux    = Data["psi"  ]
    volume  = ( zAxis[1] - zAxis[0] ) * ( rAxis[1] - rAxis[0] )
    #  -- [2-2] 流速 / 磁場 / 渦度 / ベクトルポテンシャル  --  #
    uef     = pil.pilearr( ( Data["uex"], Data["uez"], Data["uey"] ) )
    uif     = pil.pilearr( ( Data["uix"], Data["uiz"], Data["uiy"] ) )
    bfd     = pil.pilearr( ( Data["Bx" ], Data["Bz" ], Data["By" ] ) )
    wef     = crl.curlv  ( vec=uef, x1Axis=zAxis, x2Axis=rAxis, coordinate="rzt" )
    wif     = crl.curlv  ( vec=uif, x1Axis=zAxis, x2Axis=rAxis, coordinate="rzt" )
    avp     = crl.curlInv( vec=bfd, x1Axis=zAxis, x2Axis=rAxis, coordinate="rzt" )

    # --------------------------------- #
    # --- [3] ヘリシティ計算 / 返却 --- #
    # --------------------------------- #
    hlx_i   = ghl.GridHelicity( u=uif, w=wif, a=avp, b=bfd, vol=volume )
    hlx_e   = ghl.GridHelicity( u=uef, w=wef, a=avp, b=bfd, vol=volume )
    ret     = { "xAxis":zAxis, "yAxis":rAxis, "flux":flux, "HelixVi":hlx_i["HelixV"], "HelixVe":hlx_e["HelixV"], \
                "Hif":hlx_i["HelixD"][0,:,:], "Hic":hlx_i["HelixD"][1,:,:], "Him":hlx_i["HelixD"][2,:,:], \
                "Hef":hlx_e["HelixD"][0,:,:], "Hec":hlx_e["HelixD"][1,:,:], "Hem":hlx_e["HelixD"][2,:,:] }
    if ( Flag__FieldReturn ):
        ret["uex"] = uef[0,:,:] ; ret["uey"] = uef[2,:,:] ; ret["uez"] = uef[1,:,:]
        ret["uix"] = uif[0,:,:] ; ret["uiy"] = uif[2,:,:] ; ret["uiz"] = uif[1,:,:]
        ret["wex"] = wef[0,:,:] ; ret["wey"] = wef[2,:,:] ; ret["wez"] = wef[1,:,:]
        ret["wix"] = wif[0,:,:] ; ret["wiy"] = wif[2,:,:] ; ret["wiz"] = wif[1,:,:]
        ret["Bx" ] = bfd[0,:,:] ; ret["By" ] = bfd[2,:,:] ; ret["Bz" ] = bfd[1,:,:]
        ret["Avx"] = avp[0,:,:] ; ret["Avy"] = avp[2,:,:] ; ret["Avz"] = avp[1,:,:]
    return( ret )


# ------------------------------------------------------------ #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args   = rar.myRecvArgs()
    job    = args["job"   ]
    jobDir = args["jobDir"]
    kstep  = args["kstep" ]
    hlx    = FetchHelix( job=job, jobDir=jobDir, kstep=kstep, Flag__FieldReturn=True )
    print( hlx.keys() )
    # import myStyle.cMap2D as clm
    # clm.cMap2D( cMap = hlx["HelixDe"][2,:,:], Cntr=hlx["flux"] )
