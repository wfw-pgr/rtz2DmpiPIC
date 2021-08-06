import numpy              as np
import fLIB.fLIB__mpi2one as m2o
import myStyle.cMap2D     as clm

if ( __name__=='__main__' ):
    import myUtils.myRecvArgs as rar
    args  = rar.myRecvArgs()
    kstep = args["ksteps"][0]
    exe1  = m2o.mpi2one( job=args["job"], kstep=kstep )
    for i in range(18):
        FigName = "png/out{0}.png".format(i)
        print( FigName )
        print( "[mpi2cmp] {0} is under processing...".format(FigName), end="" )
        draw = clm.cMap2D( cMap=exe1[i,:,:], FigName=FigName )
        print( "\t [Done]" )
