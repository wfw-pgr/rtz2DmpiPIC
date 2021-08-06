;;; ------------------ ;;;
;;;  LoadConst         ;;;
;;; ------------------ ;;;
;;; constants.datをロードする.
Function LoadConst, job = job, config = config, GridMode = GridMode

  if not( keyword_set( config   ) ) then config   = LoadConfig()
  if not( keyword_set( job      ) ) then job      = config.job
  if not( keyword_set( GridMode ) ) then GridMode = config.GridMode
  
  jobDir   = 'job/' + job + '/' + GridMode + '/'
  Info     = { Job:Job, jobDir:jobDir }
                                ; --- parameter.dat --- ;
  line     = ''
  list     = []
  FileName = jobDir + 'parameter.dat'
  OpenR, lun, FileName, /Get_Lun
  while Not EOF( lun ) Do Begin
    readf, lun, line
    if not( ( strcmp( line, '#', 1 ) ) or ( strcmp( line, '' ) ) ) then begin
                                ; -- コメントライン( #~~ )をスキップ 
      list = [ [list], [strsplit( line, /extract )] ]
    endif
  endwhile
  free_lun, lun
  Const    = List2Structure( list = list )
  
  ret      = Create_Struct( Info, Const  )
  
  return, ret
End
