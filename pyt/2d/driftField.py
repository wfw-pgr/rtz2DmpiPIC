import numpy                 as np
import fLIB.fLIB__grad2d     as grd
import mPICsee.FetchPIC      as fpc
import myStyle.LoadConfig    as lcf
import myBasicAlgs.robustInv as riv
import myConvert.pilearr     as pil
import fLIB.fLIB__ugradv2d   as ugv
import myUtils.LoadConst     as lcn
import Filter.nxLinearFilter2D as flt

# ======================================== #
# === MHD的なドリフト場を解析          === #
# ======================================== #
def driftField( job=None, jobDir=None, kstep=None, kstep_=None, config=None, CnsFile=None, coordinate="rtz", ptype="e",
                RawData=True ):
    # ------------------------------ #
    # --- [1] 引数チェック       --- #
    # ------------------------------ #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    if ( kstep_  is None ): kstep_  = kstep - 1000
    # ------------------------------ #
    # --- [2] データ呼び出し     --- #
    # ------------------------------ #
    const   = lcn.LoadConst( InpFile=CnsFile )
    dt      = const["dt"] * float( kstep - kstep_ )
    mq      = - 1.0 if ( ptype == "e" ) else const["mr"]
    vthcv   = const["vthcv"]
    Fkeys   = [ "Bx", "By", "Bz", "Ex", "Ey", "Ez" ]
    Data1   = fpc.FetchPIC( job =job  , jobDir=jobDir, kstep  =kstep, \
                            keys=Fkeys, config=config, RawData=RawData )
    Data2   = fpc.FetchPIC( job =job  , jobDir=jobDir, kstep =kstep_, \
                            keys=Fkeys, config=config, RawData=RawData )
    xAxis   = + Data1["xAxis"]
    yAxis   = + Data1["yAxis"]
    Data1   = flt.nxLinearFilter2D( x1Axis=xAxis, x2Axis =yAxis, Data=Data1, \
                                    keys  =Fkeys, nFilter=10   , alpha=0.2   )
    Data2   = flt.nxLinearFilter2D( x1Axis=xAxis, x2Axis =yAxis, Data=Data2, \
                                    keys  =Fkeys, nFilter=10   , alpha=0.2   )

    # ------------------------------ #
    # --- [3] ドリフト速度       --- #
    # ------------------------------ #
    #  -- [3-1] 準備             --  #
    absB1   =       np.sqrt( Data1["Bx" ]**2  + Data1["By" ]**2  + Data1["Bz" ]**2    )
    absB2   =       np.sqrt( Data2["Bx" ]**2  + Data2["By" ]**2  + Data2["Bz" ]**2    )
    normB1  = pil.pilearr( ( Data1["Bx"]/absB1,  Data1["By"]/absB1,  Data1["Bz"]/absB1 ) )
    normB2  = pil.pilearr( ( Data2["Bx"]/absB2,  Data2["By"]/absB2,  Data2["Bz"]/absB2 ) )
    BsqInv1 = riv.robustInv( absB1 **2 )
    BsqInv2 = riv.robustInv( absB2 **2 )
    BcbInv1 = riv.robustInv( absB1 **3 )
    BcbInv2 = riv.robustInv( absB2 **3 )
    #  -- [3-2] ExB drift        --  #
    vExBx1  = ( Data1["Ey" ]*Data1["Bz" ] - Data1["Ez" ]*Data1["By" ] )*BsqInv1
    vExBy1  = ( Data1["Ez" ]*Data1["Bx" ] - Data1["Ex" ]*Data1["Bz" ] )*BsqInv1
    vExBz1  = ( Data1["Ex" ]*Data1["By" ] - Data1["Ey" ]*Data1["Bx" ] )*BsqInv1
    vExBx2  = ( Data2["Ey" ]*Data2["Bz" ] - Data2["Ez" ]*Data2["By" ] )*BsqInv2
    vExBy2  = ( Data2["Ez" ]*Data2["Bx" ] - Data2["Ex" ]*Data2["Bz" ] )*BsqInv2
    vExBz2  = ( Data2["Ex" ]*Data2["By" ] - Data2["Ey" ]*Data2["Bx" ] )*BsqInv2
    #  -- [3-3] gradB drift      --  #
    mqvLsq  = - 0.5 * mq * vthcv**2
    gradB1  = grd.grad2d( xg=xAxis, yg=yAxis, Data=absB1, difftype="Central__dx2" )
    gradB2  = grd.grad2d( xg=xAxis, yg=yAxis, Data=absB2, difftype="Central__dx2" )
    vgrBx1  = mqvLsq * (                             - gradB1["dfdx"]*Data1["By" ]  ) * BcbInv1
    vgrBy1  = mqvLsq * ( gradB1["dfdx"]*Data1["Bx" ] - gradB1["dfdy"]*Data1["Bz" ]  ) * BcbInv1
    vgrBz1  = mqvLsq * ( gradB1["dfdy"]*Data1["By" ]                                ) * BcbInv1
    vgrBx2  = mqvLsq * (                             - gradB2["dfdx"]*Data2["By" ]  ) * BcbInv2
    vgrBy2  = mqvLsq * ( gradB2["dfdx"]*Data2["Bx" ] - gradB2["dfdy"]*Data2["Bz" ]  ) * BcbInv2
    vgrBz2  = mqvLsq * ( gradB2["dfdy"]*Data2["By" ]                                ) * BcbInv2
    #  -- [3-4] Curvature drift  --  #
    mqvpsq  = - mq * vthcv
    bgradb1 = ugv.ugradv2d( v1=normB1, v2=normB1, zg=xAxis, xg=yAxis, coordinate=coordinate )
    bgradb2 = ugv.ugradv2d( v1=normB2, v2=normB2, zg=xAxis, xg=yAxis, coordinate=coordinate )
    vcvtx1  = mqvpsq * ( bgradb1[1,:,:]*Data1["Bz"] - bgradb1[2,:,:]*Data1["By"] ) * BsqInv1
    vcvty1  = mqvpsq * ( bgradb1[2,:,:]*Data1["Bx"] - bgradb1[0,:,:]*Data1["Bz"] ) * BsqInv1
    vcvtz1  = mqvpsq * ( bgradb1[0,:,:]*Data1["By"] - bgradb1[1,:,:]*Data1["Bx"] ) * BsqInv1
    vcvtx2  = mqvpsq * ( bgradb2[1,:,:]*Data2["Bz"] - bgradb2[2,:,:]*Data2["By"] ) * BsqInv2
    vcvty2  = mqvpsq * ( bgradb2[2,:,:]*Data2["Bx"] - bgradb2[0,:,:]*Data2["Bz"] ) * BsqInv2
    vcvtz2  = mqvpsq * ( bgradb2[0,:,:]*Data2["By"] - bgradb2[1,:,:]*Data2["Bx"] ) * BsqInv2
    #  -- [3-5] Polarization drift --  #
    #   - polarization drift -  #
    vExBx   = 0.5 * ( vExBx1 + vExBx2 )
    vExBy   = 0.5 * ( vExBy1 + vExBy2 )
    vExBz   = 0.5 * ( vExBz1 + vExBz2 )
    vgrBx   = 0.5 * ( vgrBx1 + vgrBx2 )
    vgrBy   = 0.5 * ( vgrBy1 + vgrBy2 )
    vgrBz   = 0.5 * ( vgrBz1 + vgrBz2 )
    vcvtx   = 0.5 * ( vcvtx1 + vcvtx2 )
    vcvty   = 0.5 * ( vcvty1 + vcvty2 )
    vcvtz   = 0.5 * ( vcvtz1 + vcvtz2 )
    vplrx   = ( ( vExBx1+vgrBx1+vcvtx1 ) - ( vExBx2+vgrBx2+vcvtx2 ) ) / dt
    vplry   = ( ( vExBy1+vgrBy1+vcvty1 ) - ( vExBy2+vgrBy2+vcvty2 ) ) / dt
    vplrz   = ( ( vExBz1+vgrBz1+vcvtz1 ) - ( vExBz2+vgrBz2+vcvtz2 ) ) / dt
    ret     = { "vExBx":vExBx, "vExBy":vExBy, "vExBz":vExBz, \
                "vgrBx":vgrBx, "vgrBy":vgrBy, "vgrBz":vgrBz, \
                "vcvtx":vcvtx, "vcvty":vcvty, "vcvtz":vcvtz, \
                "vplrx":vplrx, "vplry":vplry, "vplrz":vplrz  }
    return( ret )
