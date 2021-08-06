;;; ---------- ;;;
;;;  divB.pro  ;;;
;;; ---------- ;;;
;;; monopoleを可視化. ;;;
Pro divB, config = config, range = range, sfc = sfc, GridMode = GridMode, job = job

                                ; -- [1] 準備 -- ;
  if not( keyword_set( config     ) ) then config     = LoadConfig()
  if not( keyword_set( job        ) ) then job        = config.job
  if not( keyword_set( range      ) ) then range      = [0, -1, 0, -1]
  if not( keyword_set( GridMode   ) ) then GridMode   = 'RgD'
  if not( keyword_set( Coordinate ) ) then Coordinate = 'xyz'
  if    ( Coordinate eq 'xyz'       ) then keys       = [ 'Bx', 'By', 'Bz', 'psi' ]
  if    ( Coordinate eq 'rtz'       ) then keys       = [ 'Br', 'Bt', 'Bz', 'psi' ]
  const  = LoadConst( job = job, GridMode = GridMode )
  Field  = Fetch_Field( config = config, keys = keys, range = range, job = job )
  
                                ; -- [2] 計算 ( divB ) -- ;
  dx1Inv = 1.d0 / const.dz
  dx2Inv = 1.d0 / const.dr
  divB   = dblarr( const.LI, const.LJ ) * 0.d0
  if ( ( GridMode eq 'RgD' ) and ( coordinate eq 'xyz' ) ) then begin
    for j = 1, const.LJ-2 do begin
      for i = 1, const.LI-2 do begin
        divB[i, j] = + ( Field.Bx[i+1, j] - Field.Bx[i-1, j] )*dx1Inv $
                     + ( Field.By[i, j+1] - Field.By[i, j-1] )*dx2Inv
      endfor
    endfor
  endif
  RangeInfo = RangeConvert( range = range, const = const )
  print, ' Min( divB ) === ', min( divB )
  print, ' Max( divB ) === ', max( divB )
  
                                ; -- [3] 可視化 -- ;
                                ;  - [3-1] カラーマップ - ;
  if not( keyword_set( sfc ) ) then $
     My2DMap, cmp    = divB,        cnt   = Field.psi,   $
              xAxis  = Field.xAxis, yAxis = Field.yAxis, config = config
  if    ( keyword_set( sfc ) ) then $
     MySurface, data = divB, xAxis = Field.xAxis, yAxis = Field.yAxis, config = config, range = range

End
