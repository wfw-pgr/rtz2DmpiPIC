import numpy                  as np
import myStyle.plot1D         as pl1
import myStyle.LoadConfig     as lcf
import myStyle.configSettings as cfs

config   = lcf.LoadConfig()
FileName = "pwell.dat"
with open( FileName, "r" ) as f:
    Data = np.loadtxt( f, skiprows=1 )
print( Data.shape )
xAxis = Data[:,0]
yAxis = Data[:,4]

cfs.configSettings( config=config, configType="plot1D_def" )
pl1.plot1D( xAxis=xAxis, yAxis=yAxis, config=config )
