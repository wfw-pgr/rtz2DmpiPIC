import sys
import numpy              as np
import myStyle.LoadConfig as lcf

# --- 磁場 平行 / 垂直 方向に分けた座標系へと変換する --- #
def bfCoordinate( bfd=None, config=None  ):
    # ------------------------------------------- #
    # --- [1]   引数チェック                  --- #
    # ------------------------------------------- #
    if ( config  is None ): config  = lcf.LoadConfig()
    if ( bfd     is None ): sys.exit( "[bfCoordinate] ERROR!! No bfd" )
    LJ, LI       = bfd.shape[1], bfd.shape[2]

    # ------------------------------------------- #
    # --- [2] b      ( 磁場平行方向ベクトル ) --- #
    # ------------------------------------------- #
    unitbt        = np.zeros( (3,LJ,LI) )
    normb         = np.sqrt( bfd[0,:,:]**2 + bfd[1,:,:]**2 + bfd[2,:,:]**2 )
    unitbt[0,:,:] = bfd[0,:,:] / normb
    unitbt[1,:,:] = bfd[1,:,:] / normb
    unitbt[2,:,:] = bfd[2,:,:] / normb
    
    # ------------------------------------------- #
    # --- [3] n1,n2  ( 磁場垂直方向ベクトル ) --- #
    # ------------------------------------------- #
    unitbp        = np.zeros( (3,LJ,LI) )
    unitn1        = np.zeros( (3,LJ,LI) )
    unitn2        = np.zeros( (3,LJ,LI) )
    normbp        = np.sqrt( bfd[0,:,:]**2 + bfd[2,:,:]**2 )
    if ( normbp == 0.0 ): sys.exit()
    unitbp[0,:,:] = bfd[0,:,:] * normbp
    unitbp[1,:,:] = 0.0
    unitbp[2,:,:] =    bfd[2,:,:]*normbp
    unitn1[0,:,:] = unitbp[1,:,:]*unitbt[2,:,:] - unitbp[2,:,:]*unitbt[1,:,:]
    unitn1[1,:,:] = unitbp[2,:,:]*unitbt[0,:,:] - unitbp[0,:,:]*unitbt[2,:,:]
    unitn1[2,:,:] = unitbp[0,:,:]*unitbt[1,:,:] - unitbp[1,:,:]*unitbt[0,:,:]
    unitn2[0,:,:] = unitn1[1,:,:]*unitbt[2,:,:] - unitn1[2,:,:]*unitbt[1,:,:]
    unitn2[1,:,:] = unitn1[2,:,:]*unitbt[0,:,:] - unitn1[0,:,:]*unitbt[2,:,:]
    unitn2[2,:,:] = unitn1[0,:,:]*unitbt[1,:,:] - unitn1[1,:,:]*unitbt[0,:,:]

    # ------------------------------------------- #
    # --- [4] bb     ( 磁場方向 - テンソル )  --- #
    # ------------------------------------------- #
    tnsrbb        = np.zeros( (6,LJ,LI) )
    #  -   対角項 - #
    tnsrbb[0,:,:] = unitbt[0,:,:] * unitbt[0,:,:]
    tnsrbb[1,:,:] = unitbt[1,:,:] * unitbt[1,:,:]
    tnsrbb[2,:,:] = unitbt[2,:,:] * unitbt[2,:,:]
    #  - 非対角項 - #
    tnsrbb[3,:,:] = unitbt[1,:,:] * unitbt[2,:,:] # yz  ( - x )
    tnsrbb[4,:,:] = unitbt[2,:,:] * unitbt[0,:,:] # zx  ( - y )
    tnsrbb[5,:,:] = unitbt[0,:,:] * unitbt[1,:,:] # xy  ( - z )
    
    # ------------------------------------------- #
    # --- [5] I-bb   ( 磁場垂直方向テンソル ) --- #
    # ------------------------------------------- #
    tnsrIb        = np.zeros( (6,LJ,LI) )
    tnsrIb[0,:,:] = 1.0 - tnsrbb[0,:,:]
    tnsrIb[1,:,:] = 1.0 - tnsrbb[1,:,:]
    tnsrIb[2,:,:] = 1.0 - tnsrbb[2,:,:]
    tnsrIb[3,:,:] =     - tnsrbb[3,:,:]
    tnsrIb[4,:,:] =     - tnsrbb[4,:,:]
    tnsrIb[5,:,:] =     - tnsrbb[5,:,:]

    return( { "b" :unitbt, "n1"  :unitn1, "n2":unitn2, \
              "bb":tnsrbb, "Imbb":tnsrIb } )
