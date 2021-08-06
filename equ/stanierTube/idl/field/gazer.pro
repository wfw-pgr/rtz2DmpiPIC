;;; -------- ;;;
;;;  Gazer   ;;;
;;; -------- ;;;
;;; 場の量を可視化. ;;;
Pro Gazer, keys = keys, kstep = kstep, config = config, range = range $
           , job = job, sfc = sfc, psiplot = psiplot, GridMode = GridMode $
           , pdf = pdf, png = png

                                ; -- [1] 準備 -- ;
  if not( keyword_set( config    ) ) then config       = LoadConfig()
  if not( keyword_set( keys      ) ) then keys         = [ 'Jt' ]
  if not( keyword_set( GridMode  ) ) then GridMode     = config.GridMode
  if    ( keyword_set( pdf       ) ) then config.pdf   = 1
  if    ( keyword_set( png       ) ) then config.png   = 1
  if    ( n_elements( kstep ) ne 0 ) then config.kstep = kstep
  Nkeys     =  n_elements( keys )
  if    ( keyword_set( psiplot   ) ) then contline = 'psi' else contline = 'At'
  if not( KeyInclude( key = contline, array = keys ) ) then keys = [ keys, contline ]
  Field     = Fetch_Field( config = config, keys = keys, range = range, job = job,  GridMode = GridMode )
  exe       = execute( 'contline = Field.'+contline )
  
                                ; -- [1.9] 2D Surface -- ;
  if keyword_set( sfc ) then begin
    config.GraphName = keys[0]
    exe              = execute( 'cmp  = Field.' + keys[0] )
    MySurface, data = cmp, xAxis = Field.xAxis, yAxis = Field.yAxis, config = config
    return
  endif

                                ; -- [2] プロット -- ;
  for i = 0, Nkeys-1 do begin
    config.GraphName = keys[i]
    exe              = execute( 'cmp  = Field.' + keys[i] )
    My2DMap, cmp   = cmp,         cnt   = contline,   vecX = vecX,  vecY  = vecY, $
             xAxis = Field.xAxis, yAxis = Field.yAxis, config = config
  endfor


End
