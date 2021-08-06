import numpy                         as np
import mPICsee.FetchPIC              as fpc
import myStyle.LoadConfig            as lcf
import myBasicAlgs.cylindricalVolume as cvl
import mPICsee.TimeUnitConvert       as tuc
import myAnalysis.outletMask         as olm
import myAnalysis.createMask         as msk
import myAnalysis.rectangularMask    as rmk
import myUtils.LoadConst             as lcn

# ============================================ #
# ===    エネルギーインベントリー 解析     === #
# ============================================ #
def inventory( job    =None, jobDir=None, ksteps =None , CnsFile =None, \
               datFile=None, datDir=None, pngFile=None , pngDir  =None, config=None, \
               NewFile=True, MaskType="rectangular"    , Flag__integral=True, RawData=True ):
    # ------------------------------- #
    # --- [1]   引数チェック      --- #
    # ------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( ksteps  is None ): ksteps  = np.arange( config["arg_Iter"] )*config["arg_Step"]+config["arg_Init"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( datDir  is None ): datDir  = "{0}dat/".format( jobDir )
    if ( datFile is None ): datFile = datDir + "inventoryAnalysis.dat"
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    # x1lim = [ 3.0, 6.0]
    # x2lim = [-0.5,+0.5]
    x1lim = [ 6.0,12.0]
    x2lim = [-1.0,+1.0]
    # ------------------------------- #
    # --- [2]  J.E を解析         --- #
    # ------------------------------- #
    if ( NewFile ):                                # -- 新規ファイル -- #
        items = "# fstep ptime Uem Uke Uki EdJ\n"
        with open( datFile, "w"  ) as f: f.write( items )
    for i,kstep in enumerate( ksteps ):
        iEnergy = internalEnergy_calculator( job  =job, jobDir=jobDir, kstep=kstep, config=config, MaskType=MaskType, \
                                             x1lim=x1lim, x2lim=x2lim, CnsFile=CnsFile, RawData=RawData )
        eFlux   =     energyFlux_calculator( job  =job, jobDir=jobDir, kstep=kstep, config=config, MaskType=MaskType, \
                                             x1lim=x1lim, x2lim=x2lim, CnsFile=CnsFile, RawData=RawData )
        ret     = np.concatenate( [ np.array( iEnergy ), np.array(eFlux) ]  )
        with open( datFile, "ab" ) as f:
            np.savetxt( f, np.array( ret ).reshape( (1,-1) ) )
        print( ret )

        
# ============================================ #
# === エネルギーインベントリー計算ルーチン === #
# ============================================ #
def internalEnergy_calculator( job=None, jobDir=None, kstep=None, CnsFile=None, x1lim=None, x2lim=None, config=None, \
                               MaskType=None, RawData=False, Flag__DisplayMask=True ):
    # -------------------------- #
    # --- [1] データ呼び出し --- #
    # -------------------------- #
    #  -- [1-1] 読込         --  #
    fstep = float( kstep )
    ptime = tuc.TimeUnitConvert( job=job, jobDir=jobDir, kstep=kstep, unit="wpi" )
    keys  = [ "Bx","By","Bz","Ex","Ey","Ez","uex","uey","uez","uix","uiy","uiz","Jx","Jy","Jz",\
              "ne","ni","pexx","peyy","pezz","pixx","piyy","pizz","psi" ]
    config["cmp_nFilter"]    = 10
    config["cmp_LinearFilt"] = 0.2
    Data  = fpc.FetchPIC( job=job, jobDir=jobDir, kstep=kstep, keys=keys, \
                          FilterKeys=keys, RawData=RawData, config=config )
    xAxis = Data["xAxis"]
    yAxis = Data["yAxis"]
    Flux  = Data["psi"  ]
    #  -- [1-2] 準備         --  #
    if   ( MaskType == "rectangular" ):
        mask = rmk.rectangularMask( x1Axis=yAxis, x2Axis=xAxis, x1lim=x1lim, x2lim=x2lim )
    elif ( MaskType == "outlet0" ):
        mask = ( olm.outletMask( Flux=Flux ) )["outlet0"]
    elif ( MaskType == "outlet1" ):
        mask = ( olm.outletMask( Flux=Flux ) )["outlet1"]
    elif ( MaskType == "outlet2" ):
        mask = ( olm.outletMask( Flux=Flux ) )["outlet2"]
    elif ( MaskType == "LCFS"    ):
        mask = msk.createMask( Flux=Flux, Flag__LCFSMode=True )
    else:
        mask = 1.0
    volm     = cvl.cylindricalVolume( zAxis=xAxis, rAxis=yAxis )
    volMask  = volm*mask
    if ( Flag__DisplayMask ):
        FigName = "mask_{0:08}.png".format( kstep )
        import myStyle.cMap2D as clm
        clm.cMap2D( xAxis=yAxis, yAxis=xAxis, cMap=mask, Cntr=Flux, FigName=FigName, config=config )
    # ----------------------------- #
    # --- [2]  Internal Energy  --- #
    # ----------------------------- #
    const    = lcn.LoadConst( InpFile=CnsFile )
    EM       =   0.5 * ( Data["Bx"]**2 + Data["By"]**2 + Data["Bz"]**2 + \
                         Data["Ex"]**2 + Data["Ey"]**2 + Data["Ez"]**2 )
    Ke       = ( 0.5 * const["rme"] * Data["ne"] * ( Data["uex"]**2 + Data["uey"]**2 + Data["uez"]**2 ) +  \
                 0.5 * ( Data["pexx"] + Data["peyy"] + Data["pezz"] ) ) * const["wpewce"]**2
    Ki       = ( 0.5 * const["rmi"] * Data["ni"] * ( Data["uix"]**2 + Data["uiy"]**2 + Data["uiz"]**2 ) + \
                 0.5 * ( Data["pixx"] + Data["piyy"] + Data["pizz"] ) ) * const["wpewce"]**2
    EJ       = Data["Ex"]*Data["Jx"] + Data["Ey"]*Data["Jy"] + Data["Ez"]*Data["Jz"]
    Uem      = np.sum( volMask * EM )
    Uke      = np.sum( volMask * Ke )
    Uki      = np.sum( volMask * Ki )
    Uej      = np.sum( volMask * EJ )
    # ----------------------------- #
    # --- [3]     返却          --- #
    # ----------------------------- #
    return( [ fstep, ptime, Uem, Uke, Uki, Uej ] )


def energyFlux_calculator( job=None, jobDir=None, kstep=None, config=None, \
                           MaskType=None, RawData=False, CnsFile=None, x1lim=None, x2lim=None ):
    # -------------------------- #
    # --- [1] データ呼び出し --- #
    # -------------------------- #
    #  -- [1-1] 読込         --  #
    fstep  = float( kstep )
    ptime  = tuc.TimeUnitConvert( job=job, jobDir=jobDir, kstep=kstep, unit="wpi" )
    keys   = [ "Bx","By","Bz","Ex","Ey","Ez","uex","uey","uez","uix","uiy","uiz","psi",\
               "ne","ni","pexx","peyy","pezz","pixx","piyy","pizz","pexy","peyz","pexz","pixy","piyz","pixz" ]
    config["cmp_nFilter"]    = 10
    config["cmp_LinearFilt"] = 0.2
    Data   = fpc.FetchPIC( job=job, jobDir=jobDir, kstep=kstep, keys =keys, \
                           FilterKeys=keys, RawData=RawData, config  =config )
    x1Axis = Data["yAxis"]
    x2Axis = Data["xAxis"]
    # ---------------------------------------- #
    # --- [2] 準備                         --- #
    # ---------------------------------------- #
    const  = lcn.LoadConst( InpFile=CnsFile )
    dr,dz  = x1Axis[1]-x1Axis[0], x2Axis[1]-x2Axis[0]
    dij    = 2
    jp1    = np.argmin( np.abs( x1Axis - x1lim[0] ) )
    jp2    = np.argmin( np.abs( x1Axis - x1lim[1] ) )
    ip1    = np.argmin( np.abs( x2Axis - x2lim[0] ) )
    ip2    = np.argmin( np.abs( x2Axis - x2lim[1] ) )
    rp     = x1Axis[(jp1):(jp2+1)]
    dSnz1  = 2.0 * np.pi * rp * dr * ( -1.0 )
    dSnz2  = 2.0 * np.pi * rp * dr * ( +1.0 )
    dSnr1  = 2.0 * np.pi * x1Axis[jp1] * dz * ( -1.0 )
    dSnr2  = 2.0 * np.pi * x1Axis[jp2] * dz * ( +1.0 )
    # ---------------------------------------- #
    # --- [3] Poyinting Flux               --- #
    # ---------------------------------------- #
    PoyFx  = Data["Ey"]*Data["Bz"] - Data["Ez"]*Data["By"]
    PoyFz  = Data["Ex"]*Data["By"] - Data["Ey"]*Data["Bx"]
    dPoyFx = np.sum( np.mean( PoyFx[(jp1-dij):(jp1+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr1 )
    uPoyFx = np.sum( np.mean( PoyFx[(jp2-dij):(jp2+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr2 )
    lPoyFz = np.sum( np.mean( PoyFz[(jp1):(jp2+1),(ip1-dij):(ip1+dij+1)], axis=1 ) * dSnz1 )
    rPoyFz = np.sum( np.mean( PoyFz[(jp1):(jp2+1),(ip2-dij):(ip2+dij+1)], axis=1 ) * dSnz2 )
    delPoyF= dPoyFx + uPoyFx + lPoyFz + rPoyFz
    # ---------------------------------------- #
    # --- [4] Particle Energy Flux (elec.) --- #
    # ---------------------------------------- #
    Ket     = 0.5 * const["wpewce"]**2 * const["rme"] * Data["ne"] * ( Data["uex"]**2 + Data["uey"]**2 + Data["uez"]**2 )
    pet     = 0.5 * const["wpewce"]**2 * ( Data["pexx"] + Data["peyy"] + Data["pezz"] )
    #  -- [4-1] 運動エネルギー移流     -- #
    Kex     = Ket * Data["uex"]
    Kez     = Ket * Data["uez"]
    dKex    = np.sum( np.mean(    Kex[(jp1-dij):(jp1+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr1 )
    uKex    = np.sum( np.mean(    Kex[(jp2-dij):(jp2+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr2 )
    lKez    = np.sum( np.mean(    Kez[(jp1):(jp2+1),(ip1-dij):(ip1+dij+1)], axis=1 ) * dSnz1 )
    rKez    = np.sum( np.mean(    Kez[(jp1):(jp2+1),(ip2-dij):(ip2+dij+1)], axis=1 ) * dSnz2 )
    #  -- [4-2] 熱エネルギー移流       -- #
    peux    = pet * Data["uex"]
    peuz    = pet * Data["uez"]
    dpeux   = np.sum( np.mean(   peux[(jp1-dij):(jp1+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr1 )
    upeux   = np.sum( np.mean(   peux[(jp2-dij):(jp2+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr2 )
    lpeuz   = np.sum( np.mean(   peuz[(jp1):(jp2+1),(ip1-dij):(ip1+dij+1)], axis=1 ) * dSnz1 )
    rpeuz   = np.sum( np.mean(   peuz[(jp1):(jp2+1),(ip2-dij):(ip2+dij+1)], axis=1 ) * dSnz2 )
    #  -- [4-3] 対角項熱エネルギー移流 -- #
    dgpeux  = Data["pexx"] * Data["uex"] * const["wpewce"]**2
    dgpeuz  = Data["pezz"] * Data["uez"] * const["wpewce"]**2
    ddgpeux = np.sum( np.mean( dgpeux[(jp1-dij):(jp1+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr1 )
    udgpeux = np.sum( np.mean( dgpeux[(jp2-dij):(jp2+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr2 )
    ldgpeuz = np.sum( np.mean( dgpeuz[(jp1):(jp2+1),(ip1-dij):(ip1+dij+1)], axis=1 ) * dSnz1 )
    rdgpeuz = np.sum( np.mean( dgpeuz[(jp1):(jp2+1),(ip2-dij):(ip2+dij+1)], axis=1 ) * dSnz2 )
    #  -- [4-4] 対角項熱エネルギー移流 -- #
    gvteux  = ( Data["pexy"] * Data["uey"] + Data["pexz"] * Data["uez"] ) * const["wpewce"]**2
    gvteuz  = ( Data["pexz"] * Data["uex"] + Data["peyz"] * Data["uey"] ) * const["wpewce"]**2
    dgvteux = np.sum( np.mean( gvteux[(jp1-dij):(jp1+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr1 )
    ugvteux = np.sum( np.mean( gvteux[(jp2-dij):(jp2+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr2 )
    lgvteuz = np.sum( np.mean( gvteuz[(jp1):(jp2+1),(ip1-dij):(ip1+dij+1)], axis=1 ) * dSnz1 )
    rgvteuz = np.sum( np.mean( gvteuz[(jp1):(jp2+1),(ip2-dij):(ip2+dij+1)], axis=1 ) * dSnz2 )
    #  -- [4-5] 熱フラックス           -- #
    qeFx    = 0.0
    qeFz    = 0.0
    dqeFx   = 0.0 * qeFx
    uqeFx   = 0.0 * qeFx
    lqeFz   = 0.0 * qeFz
    rqeFz   = 0.0 * qeFz
    delKe   =    dKex +    uKex +    lKez +    rKez
    delpeu  =   dpeux +   upeux +   lpeuz +   rpeuz
    deldgpeu= ddgpeux + udgpeux + ldgpeuz + rdgpeuz 
    delgvteu= dgvteux + ugvteux + lgvteuz + rgvteuz 
    delqeF  =   dqeFx +   uqeFx +   lqeFz +   rqeFz
    detot   = dKex +  dpeux +  ddgpeux +  dgvteux +  dqeFx
    uetot   = uKex +  upeux +  udgpeux +  ugvteux +  uqeFx
    letot   = lKez +  lpeuz +  ldgpeuz +  lgvteuz +  lqeFz
    retot   = rKez +  rpeuz +  rdgpeuz +  rgvteuz +  rqeFz
    eFtotal = delKe+ delpeu + deldgpeu + delgvteu + delqeF

    # dqFluxx = np.sum( np.mean( qFluxx[(jp1-dij):(jp1+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr1 )
    # uqFluxx = np.sum( np.mean( qFluxx[(jp2-dij):(jp2+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr2 )
    # lqFluxz = np.sum( np.mean( qFluxz[(jp1):(jp2+1),(ip1-dij):(ip1+dij+1)], axis=1 ) * dSnz1 )
    # rqFluxz = np.sum( np.mean( qFluxz[(jp1):(jp2+1),(ip2-dij):(ip2+dij+1)], axis=1 ) * dSnz2 )    
    # Key    = Ket * Data["uey"]
    # peuy   = pet * Data["uey"]
    # qFluxy = 0.0
    # gvteuy = Data["pexy"] * Data["uex"] + Data["peyz"] * Data["uez"]
    # dgpeuy = Data["peyy"] * Data["uey"]
    # ---------------------------------------- #
    # --- [5] Particle Energy Flux ( ion ) --- #
    # ---------------------------------------- #
    Kit     = 0.5 * const["wpewce"]**2 * const["rmi"] * Data["ni"] * ( Data["uix"]**2 + Data["uiy"]**2 + Data["uiz"]**2 )
    pit     = 0.5 * const["wpewce"]**2 * ( Data["pixx"] + Data["piyy"] + Data["pizz"] )
    #  -- [5-1] 運動エネルギー移流     -- #
    Kix     = Kit * Data["uix"]
    Kiz     = Kit * Data["uiz"]
    dKix    = np.sum( np.mean(    Kix[(jp1-dij):(jp1+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr1 )
    uKix    = np.sum( np.mean(    Kix[(jp2-dij):(jp2+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr2 )
    lKiz    = np.sum( np.mean(    Kiz[(jp1):(jp2+1),(ip1-dij):(ip1+dij+1)], axis=1 ) * dSnz1 )
    rKiz    = np.sum( np.mean(    Kiz[(jp1):(jp2+1),(ip2-dij):(ip2+dij+1)], axis=1 ) * dSnz2 )
    #  -- [5-2] 熱エネルギー移流       -- #
    piux    = pit * Data["uix"]
    piuz    = pit * Data["uiz"]
    dpiux   = np.sum( np.mean(   piux[(jp1-dij):(jp1+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr1 )
    upiux   = np.sum( np.mean(   piux[(jp2-dij):(jp2+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr2 )
    lpiuz   = np.sum( np.mean(   piuz[(jp1):(jp2+1),(ip1-dij):(ip1+dij+1)], axis=1 ) * dSnz1 )
    rpiuz   = np.sum( np.mean(   piuz[(jp1):(jp2+1),(ip2-dij):(ip2+dij+1)], axis=1 ) * dSnz2 )
    #  -- [5-3] 対角項熱エネルギー移流 -- #
    dgpiux  = ( Data["pixx"] * Data["uix"] ) * const["wpewce"]**2
    dgpiuz  = ( Data["pizz"] * Data["uiz"] ) * const["wpewce"]**2
    ddgpiux = np.sum( np.mean( dgpiux[(jp1-dij):(jp1+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr1 )
    udgpiux = np.sum( np.mean( dgpiux[(jp2-dij):(jp2+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr2 )
    ldgpiuz = np.sum( np.mean( dgpiuz[(jp1):(jp2+1),(ip1-dij):(ip1+dij+1)], axis=1 ) * dSnz1 )
    rdgpiuz = np.sum( np.mean( dgpiuz[(jp1):(jp2+1),(ip2-dij):(ip2+dij+1)], axis=1 ) * dSnz2 )
    #  -- [5-4] 対角項熱エネルギー移流 -- #
    gvtiux  = ( Data["pixy"] * Data["uiy"] + Data["pixz"] * Data["uiz"] ) * const["wpewce"]**2
    gvtiuz  = ( Data["pixz"] * Data["uix"] + Data["piyz"] * Data["uiy"] ) * const["wpewce"]**2
    dgvtiux = np.sum( np.mean( gvtiux[(jp1-dij):(jp1+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr1 )
    ugvtiux = np.sum( np.mean( gvtiux[(jp2-dij):(jp2+dij+1),(ip1):(ip2+1)], axis=0 ) * dSnr2 )
    lgvtiuz = np.sum( np.mean( gvtiuz[(jp1):(jp2+1),(ip1-dij):(ip1+dij+1)], axis=1 ) * dSnz1 )
    rgvtiuz = np.sum( np.mean( gvtiuz[(jp1):(jp2+1),(ip2-dij):(ip2+dij+1)], axis=1 ) * dSnz2 )
    #  -- [5-5] 熱フラックス           -- #
    qiFx    = 0.0
    qiFz    = 0.0
    dqiFx   = 0.0 * qiFx
    uqiFx   = 0.0 * qiFx
    lqiFz   = 0.0 * qiFz
    rqiFz   = 0.0 * qiFz
    delKi   =    dKix +    uKix +    lKiz +    rKiz
    delpiu  =   dpiux +   upiux +   lpiuz +   rpiuz
    deldgpiu= ddgpiux + udgpiux + ldgpiuz + rdgpiuz 
    delgvtiu= dgvtiux + ugvtiux + lgvtiuz + rgvtiuz 
    delqiF  =   dqiFx +   uqiFx +   lqiFz +   rqiFz
    ditot   = dKix +  dpiux +  ddgpiux +  dgvtiux +  dqiFx
    uitot   = uKix +  upiux +  udgpiux +  ugvtiux +  uqiFx
    litot   = lKiz +  lpiuz +  ldgpiuz +  lgvtiuz +  lqiFz
    ritot   = rKiz +  rpiuz +  rdgpiuz +  rgvtiuz +  rqiFz
    iFtotal = delKi+ delpiu + deldgpiu + delgvtiu + delqiF
    # ----------------------------- #
    # --- [6]     返却          --- #
    # ----------------------------- #
    print( " === Energy Inventory === " )
    print( "[Poynting]" )
    print( "\t{0:>+.6f}  {1:>+.6f}  {2:>+.6f}  {3:>+.6f}  {4:>+.6f}".format(  dPoyFx,  uPoyFx,  lPoyFz,  rPoyFz,  delPoyF ) )
    print( "[electron]" )
    print( "\t K  \t{0:>+.6f}  {1:>+.6f}  {2:>+.6f}  {3:>+.6f}  {4:>+.6f}".format(    dKex,    uKex,    lKez,    rKez,    delKe ) )
    print( "\t pu \t{0:>+.6f}  {1:>+.6f}  {2:>+.6f}  {3:>+.6f}  {4:>+.6f}".format(   dpeux,   upeux,   lpeuz,   rpeuz,   delpeu ) )
    print( "\t dp \t{0:>+.6f}  {1:>+.6f}  {2:>+.6f}  {3:>+.6f}  {4:>+.6f}".format( ddgpeux, udgpeux, ldgpeuz, rdgpeuz, deldgpeu ) )
    print( "\t gv \t{0:>+.6f}  {1:>+.6f}  {2:>+.6f}  {3:>+.6f}  {4:>+.6f}".format( dgvteux, ugvteux, lgvteuz, rgvteuz, delgvteu ) )
    print( "\t qF \t{0:>+.6f}  {1:>+.6f}  {2:>+.6f}  {3:>+.6f}  {4:>+.6f}".format(   dqeFx,   uqeFx,   lqeFz,   rqeFz,   delqeF ) )
    print( "-----------------------------------------------------------------------------" )
    print( "\t qF \t{0:>+.6f}  {1:>+.6f}  {2:>+.6f}  {3:>+.6f}  {4:>+.6f}".format(   detot,   uetot,   letot,   retot,  eFtotal ) )
    print( "[ion]" )
    print( "\t K  \t{0:>+.6f}  {1:>+.6f}  {2:>+.6f}  {3:>+.6f}  {4:>+.6f}".format(    dKix,    uKix,    lKiz,    rKiz,    delKi ) )
    print( "\t pu \t{0:>+.6f}  {1:>+.6f}  {2:>+.6f}  {3:>+.6f}  {4:>+.6f}".format(   dpiux,   upiux,   lpiuz,   rpiuz,   delpiu ) )
    print( "\t dp \t{0:>+.6f}  {1:>+.6f}  {2:>+.6f}  {3:>+.6f}  {4:>+.6f}".format( ddgpiux, udgpiux, ldgpiuz, rdgpiuz, deldgpiu ) )
    print( "\t gv \t{0:>+.6f}  {1:>+.6f}  {2:>+.6f}  {3:>+.6f}  {4:>+.6f}".format( dgvtiux, ugvtiux, lgvtiuz, rgvtiuz, delgvtiu ) )
    print( "\t qF \t{0:>+.6f}  {1:>+.6f}  {2:>+.6f}  {3:>+.6f}  {4:>+.6f}".format(   dqiFx,   uqiFx,   lqiFz,   rqiFz,   delqiF ) )
    print( "-----------------------------------------------------------------------------" )
    print( "\t qF \t{0:>+.6f}  {1:>+.6f}  {2:>+.6f}  {3:>+.6f}  {4:>+.6f}".format(   ditot,   uitot,   litot,   ritot,  iFtotal ) )
    return( [ fstep, ptime, delPoyF, delKe, delpeu, deldgpeu, delgvteu, delqeF ] )

# ======================================== #
# ===           テスト実行             === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args     = rar.myRecvArgs()
    OutDir   = "./png/inventory/"
    inventory( job   =args["job"]   , jobDir=args["jobDir"], \
                 ksteps=args["ksteps"]  )
