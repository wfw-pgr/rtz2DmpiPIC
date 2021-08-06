import sys
import numpy              as np
import myStyle.LoadConfig as lcf
import myConvert.pilearr  as pil
import myStyle.myBarPlot  as mbp
import myStyle.pickColor  as pcl

# --- 各エネルギーの構成 Barplot  --- #
def EnergyComponents( job1   =None, ptcFile1=None, totFile1=None, fldFile1=None, \
                      job2   =None, ptcFile2=None, totFile2=None, fldFile2=None, \
                      OutFile=None, OutDir  =None, config  =None  ):
    # --- [1] 引数チェック --- #
    if ( config   is     None ): config   = lcf.LoadConfig()
    if ( job1     is     None ): sys.exit(" [ERROR] -@EnergyBeforeAfter.py- job1 ???")
    if ( job2     is     None ): sys.exit(" [ERROR] -@EnergyBeforeAfter.py- job2 ???")
    if ( ptcFile1 is     None ): ptcFile1 = "../../../job/{0}/dat/energy_kinetic.dat".format( job1 )
    if ( fldFile1 is     None ): fldFile1 = "../../../job/{0}/dat/energy_field.dat".format( job1 )
    if ( totFile1 is     None ): totFile1 = "../../../job/{0}/dat/energy.dat".format( job1 )
    if ( ptcFile2 is     None ): ptcFile2 = "../../../job/{0}/dat/energy_kinetic.dat".format( job2 )
    if ( fldFile2 is     None ): fldFile2 = "../../../job/{0}/dat/energy_field.dat".format( job2 )
    if ( totFile2 is     None ): totFile2 = "../../../job/{0}/dat/energy.dat".format( job2 )
    if ( OutFile  is     None ): OutFile  = "EnergyComponents.png"
    if ( OutDir   is not None ): OutFile  = OutDir + OutFile
    
    # --- [2] エネルギー読込 --- #
    #   - job1   - #
    Eknt1    = ( np.loadtxt( ptcFile1 ) )[:,:]
    Efld1    = ( np.loadtxt( fldFile1 ) )[:,:]
    closeIdx = np.argmin( np.abs( Efld1[:,1] - 2000.0 ) )
    Uall1    = pil.pilearr( ( Efld1[closeIdx,2:5], Efld1[closeIdx,6:9], \
                              Eknt1[closeIdx,2:5], Eknt1[closeIdx,6:9] ), \
                            NoNewAxis=True, axis=0 )
    #   - job2   - #
    Eknt2    = ( np.loadtxt( ptcFile2 ) )[:,:]
    Efld2    = ( np.loadtxt( fldFile2 ) )[:,:]
    closeIdx = np.argmin( np.abs( Efld2[:,1] - 2000.0 ) )
    Uall2    = pil.pilearr( ( Efld2[closeIdx,2:5], Efld2[closeIdx,6:9], \
                              Eknt2[closeIdx,2:5], Eknt2[closeIdx,6:9] ), \
                            NoNewAxis=True, axis=0 )
    #   - job0   - #
    Uall0    = pil.pilearr( ( Efld2[0,2:5], Efld2[0,6:9], \
                              Eknt2[0,2:5], Eknt2[0,6:9] ), \
                            NoNewAxis=True, axis=0 )
    #   - 連結   - #
    Uall     = pil.pilearr( (Uall0, Uall1,Uall2) )
    #   - 規格化 - #
    for i in range( Uall.shape[0] ):
        Uall[i,:] = Uall[i,:] / np.sum( Uall[i,:] )

    # --- [3] Barplot で プロット --- #
    #  -- コンフィグ -- #
    config["FigSize"]         = (4.5,3)
    config["plt_Position"]    = [0.02,0.25,0.75,0.98]
    config["plt_LegBoxPosit"] = [1.20,0.9,1.30,1.05]
    config["plt_LegFontSize"] = 15
    config["plt_LegNColumn"]  = 1
    config["chsize1"]         = 50
    config["chsize2"]         = 50
    
    #  -- プロット   -- #
    xlabel = ["Initial", "Case I", "Case O"]
    labels = [ "$E_r$" , "$E_t$" , "$E_z$" , "$B_r$" , "$B_t$" , "$B_z$" ,
               "$v_{er}$", "$v_{et}$", "$v_{ez}$", "$v_{ir}$", "$v_{it}$", "$v_{iz}$" ]
    colors = pcl.pickColor( pallete="EBvevi", nColor=12 )
    fig    = mbp.myBarPlot( Data    =Uall  , \
                            StackBar= True , \
                            xLabel  =xlabel, \
                            labels  =labels, \
                            colors  =colors, \
                            FigName =OutFile, \
                            xLabelRotation=35, \
                            config  =config )
    print( "{0} is processed...".format(OutFile) )


# --------------------------------------- #
if ( __name__=="__main__" ):
    import myUtils.myRecvArgs as rar
    args      = rar.myRecvArgs()
    job1      = "ctrI02_"
    job2      = "ctrO02_"
    OutDir    = "./png/"
    Output    = EnergyComponents( job1=job1, job2=job2, OutDir=OutDir )
