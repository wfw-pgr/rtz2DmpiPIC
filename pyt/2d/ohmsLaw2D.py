import myStyle.cMap2D         as clm
import mPICsee.FetchPIC       as fpc
import myStyle.LoadConfig     as lcf
import myUtils.LoadConst      as lcn
import fLIB.fLIB__grad2d      as grd
import myBasicAlgs.robustInv  as inv
import myBasicAlgs.makegrid   as mkg
import Filter.LinearFilter2D  as flt
import myStyle.configSettings as cfs

# ======================================== #
# ===  Ohm's Law を 評価               === #
# ======================================== #
def ohmsLaw2D( job  =None, jobDir=None, kstep  =None, OutFile=None, OutDir=None, \
               keys =None, config=None, CnsFile=None ):
    # ------------------------ #
    # --- [1] 引数チェック --- #
    # ------------------------ #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( job     is None ): job     = config["pic_job"]
    if ( kstep   is None ): kstep   = config["pic_kstep"]
    if ( jobDir  is None ): jobDir  = "{0}{1}/".format( config["pic_jobDir"], job )
    if ( OutDir  is None ): OutDir  = "{0}png/".format( jobDir )
    if ( OutFile is None ): OutFile = "OhmsLaw2D.png"
    if ( CnsFile is None ): CnsFile = "{0}dat/constants.dat".format( jobDir )
    OutFile = OutDir + OutFile
    
    # -------------------------- #
    # --- [2] データ呼び出し --- #
    # -------------------------- #
    const  = lcn.LoadConst( InpFile=CnsFile )
    keys   = ["Ey", "uix", "uiy", "uiz" , "uex" , "uey", "uez" , "psi" , \
              "Bx", "By" , "Bz" , "pexy", "peyz", "ne" , "pixy", "piyz", "ni" ]
    Data   = fpc.FetchPIC( job =job , jobDir=jobDir, kstep =kstep, keys=keys, config=config )
    xAxis  = + Data["xAxis"]
    yAxis  = + Data["yAxis"]
    flux   = + Data["psi"  ]
    yInv   = inv.robustInv( yAxis )
    rIg    = ( mkg.makegrid( x1=yInv, x2=xAxis, indexing='xy' ) )["x1g"]
    Ey     = Data["Ey"]
    qe, qi = -1.0, +1.0

    # ------------------------------ #
    # --- [3] Ohm's Law ( 電子 ) --- #
    # ------------------------------ #
    # --    -uexB     -- #
    muexB  = - ( Data["uez"]*Data["Bx"] - Data["uex"]*Data["Bz"] )
    # -- (ue.grad) ue -- #
    grduet = grd.grad2d( Data=Data["uey"], xg=xAxis, yg=yAxis  )
    ugrdue = ( Data["uez"]*grduet["dfdx"] + Data["uex"]*grduet["dfdy"] + Data["uey"]* Data["uex" ] * rIg ) * ( qe * const["rme"] )
    # --    div pe    -- #
    qneInv = inv.robustInv( Data["ne"], Flag__positive=True  ) / qe
    divPrt = ( grd.grad2d( Data=Data["pexy"], xg=xAxis, yg=yAxis ) )["dfdy"]
    divPzt = ( grd.grad2d( Data=Data["peyz"], xg=xAxis, yg=yAxis ) )["dfdx"]
    divPe  = ( divPrt + divPzt )*qneInv
    # ------------------------------ #
    # --- [4] Ohm's Law (イオン) --- #
    # ------------------------------ #
    # -- uixB -- #
    muixB  = - ( Data["uiz"]*Data["Bx"] - Data["uix"]*Data["Bz"] )
    # -- (ui.grad) ui -- #
    grduit = grd.grad2d( Data=Data["uiy"], xg=xAxis, yg=yAxis  )
    ugrdui = ( Data["uiz"]*grduit["dfdx"] + Data["uix"]*grduit["dfdy"] + Data["uiy"]* Data["uix" ] * rIg ) * ( qi * const["rmi"] )
    # -- div pi -- #
    qniInv = inv.robustInv( Data["ni"], Flag__positive=True ) / qi
    divPrt = ( grd.grad2d( Data=Data["pixy"], xg=xAxis, yg=yAxis ) )["dfdy"]
    divPzt = ( grd.grad2d( Data=Data["piyz"], xg=xAxis, yg=yAxis ) )["dfdx"]
    divPi  = ( divPrt + divPzt )*qniInv

    # ------------------------- #
    # --- [5] プロット      --- #
    # ------------------------- #
    cfs.configSettings( config=config, configType="cMap2D_thesis" )
    nFilt  = 5
    config["cmp_LinearFilt"] = 0.5
    config["cmp_AutoLevel"]  = False
    config["cmp_MaxMin"]     = [-0.1,+0.1]
    config["cmp_ColorTable"] = "bwr"
    cMaps  = { "Ey":Ey, "uexB":muexB, "ugrdue":ugrdue, "divPe":divPe, "uixB":muixB, "ugrdui":ugrdui, "divPi":divPi }
    pKeys  = [ "Ey", "uexB", "ugrdue", "divPe", "uixB", "ugrdui", "divPi" ]
    for key in pKeys:
        cMap    = cMaps[ key ]
        FigName = OutFile.replace( ".png", "_{0}_{1:08}.png".format( key, kstep ) )
        for i in range( nFilt ):
            cMap = flt.LinearFilter2D( Data=cMap   , alpha=config["cmp_LinearFilt"], \
                                       x1Axis=xAxis, x2Axis=yAxis, coordinate="rtz" )
        clm.cMap2D( xAxis  =yAxis  , yAxis =xAxis , cMap=cMap, Cntr=flux, \
                    FigName=FigName, config=config )
        print( "[ohmsLaw2D] {0} is outputed...".format( FigName ) )

        
# ======================================== #
# ===  実行部分                        === #
# ======================================== #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args    = rar.myRecvArgs()
    OutDir  = "./png/ohmsLaw2D/"
    for kstep in args["ksteps"]:
        ohmsLaw2D( job=args["job"], jobDir=args["jobDir"], kstep=kstep, OutDir=OutDir )
