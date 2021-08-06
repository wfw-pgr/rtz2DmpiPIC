;;; ------------- ;;;
;;;  Fetch_Field  ;;;
;;; ------------- ;;;
;;; データ取得用ルーチン ;;;
Function Fetch_Field, keys = keys, range = range, config = config, job = job, GridMode = GridMode

                                ;  --  [1]  準備  -- ;
                                ;  --  [1-1] データ読込  -- ;
  if not( keyword_set( config   ) ) then config   = LoadConfig()
  if not( keyword_set( range    ) ) then range    = config.range
  if not( keyword_set( keys     ) ) then keys     = [ 'psi' ]
  if not( keyword_set( job      ) ) then job      = config.job
  if not( keyword_set( GridMode ) ) then GridMode = config.GridMode
  const     = LoadConst( job = job, GridMode = GridMode )
  LI        = const.LI
  LJ        = const.LJ
                                ;  --  [1-2] プロット範囲  -- ;
  RangeInfo = RangeConvert( range = range, const = const )
  xNum      = RangeInfo.xNum
  yNum      = RangeInfo.yNum
  xRange    = RangeInfo.xRange
  yRange    = RangeInfo.yRange
  
                                ;  --  [2] データ読み込み -- ;
                                ;  --  [2-1] 場の量 ファイル読込  -- ;
  filename  = const.JobDir + '/frp.bin'
  frp       = MyLoadByn( fname = filename, N = 1, size = [LI, LJ, 3], type = 'double' )
  filename  = const.JobDir + '/Bfd.bin'
  Bfd       = MyLoadByn( fname = filename, N = 1, size = [3, LI, LJ], type = 'double' )
  filename  = const.JobDir + '/Jcr.bin'
  Jcr       = MyLoadByn( fname = filename, N = 1, size = [3, LI, LJ], type = 'double' )
  filename  = const.JobDir + '/uvc.bin'
  uvc       = MyLoadByn( fname = filename, N = 1, size = [3, LI, LJ], type = 'double' )
                                ;  --  [2-2] 必要変数の代入   -- ;
                                ;   -  [2-2-1] RTZ 円筒座標系  - ;
  if ( KeyInclude( array = keys, key = 'Br' ) ) then Br  = reform( Bfd[ 0, xRange, yRange], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'Bt' ) ) then Bt  = reform( Bfd[ 1, xRange, yRange], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'Bz' ) ) then Bz  = reform( Bfd[ 2, xRange, yRange], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'Jr' ) ) then Jr  = reform( Jcr[ 0, xRange, yRange], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'Jt' ) ) then Jt  = reform( Jcr[ 1, xRange, yRange], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'Jz' ) ) then Jz  = reform( Jcr[ 2, xRange, yRange], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'ur' ) ) then ur  = reform( uvc[ 0, xRange, yRange], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'ut' ) ) then ut  = reform( uvc[ 1, xRange, yRange], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'uz' ) ) then uz  = reform( uvc[ 2, xRange, yRange], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'At' ) ) then At  = reform( frp[ xRange, yRange, 0], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'psi') ) then psi = reform( frp[ xRange, yRange, 0], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'rho') ) then rho = reform( frp[ xRange, yRange, 1], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'prs') ) then prs = reform( frp[ xRange, yRange, 2], xNum, yNum )
                                ;   -  [2-2-2] xyz 直交座標系  - ;
  if ( KeyInclude( array = keys, key = 'Bx' ) ) then Bx  = reform( Bfd[ 0, xRange, yRange], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'By' ) ) then By  = reform( Bfd[ 1, xRange, yRange], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'Jx' ) ) then Jx  = reform( Jcr[ 0, xRange, yRange], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'Jy' ) ) then Jy  = reform( Jcr[ 1, xRange, yRange], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'ux' ) ) then ux  = reform( uvc[ 0, xRange, yRange], xNum, yNum )
  if ( KeyInclude( array = keys, key = 'uy' ) ) then uy  = reform( uvc[ 1, xRange, yRange], xNum, yNum )

                                ;  --  [2-3] 座標軸 -- ;
  xAxis     = ( const.dz * ( dindgen( const.LI ) - 1.d0 ) )[xRange]
  yAxis     = ( const.dr * ( dindgen( const.LJ ) - 1.d0 ) )[yRange]
                                ;  --  [2-4] psi   -- ;
  if ( KeyInclude( array = keys, key = 'psi') ) then $
     psi = psi* ( replicate( 1.0, n_elements(xAxis) ) # yAxis )
  
                                ;  --  [3] 返却用構造体の生成  --  ;
  retcmd    = ' const:const, xAxis:xAxis, yAxis:yAxis'
  for i = 0, n_elements( keys )-1 do $
     retcmd = retcmd + ', ' + keys[i] + ':' + keys[i]
  ret       = {}
  retcmd    = 'ret = {' + retcmd + ' }'
  print, retcmd
  exe       = execute( retcmd )
  
  return, ret  
End
