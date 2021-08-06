;;; -------- ;;;
;;;  Gazer   ;;;
;;; -------- ;;;
;;; 場の量を可視化. ;;;
Pro Gazer, keys = keys, kstep = kstep, config = config, range = range, job = job, dual = dual

                                ; -- [1] 準備 -- ;
  if not( keyword_set( config ) ) then config       = LoadConfig()
  if not( keyword_set( keys   ) ) then keys         = [ 'Jt' ]
  if ( n_elements( kstep ) ne 0 ) then config.kstep = kstep
  keys      = [ keys, 'psi' ]
  Field     = Fetch_Field( config = config, keys = keys, range = range, job = job )
  
                                ; -- [2] プロット -- ;
  for i = 0, n_elements( keys )-2 do begin
    config.GraphName = keys[i]
    exe              = execute( 'cmp  = Field.' + keys[i] )
    Mysurface, data   = cmp, $
               xAxis = Field.zAxis, yAxis = Field.rAxis, config = config
  endfor
  
End
