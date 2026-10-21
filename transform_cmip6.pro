
;transform
;
pro transform_cmip6

  start_time=systime(1)

  infile = 'E:\RS_Data\Trendy_2026\RA\Output_Annual1\RA_ISAM_S3_2001_2013_M2.dat'

  column = 720
  row = 360
  band = 156

  data = read_binary(infile,data_dims = [column,row,band],data_type = 5)

  result = make_array(column, row, band, type = 5)

  result[0:359,*,*] = data[360:719,0:359,*]
  result[360:719,*,*] = data[0:359,0:359,*]

  out_file = strcompress('E:\RS_Data\Trendy_2026\RA\Output_Annual1\'+'RA_ISAM_S3_2001_2013_M'+'.dat',/remove_all)
  openw,out_lun,out_file,/get_lun
  writeu,out_lun,result
  free_lun,out_lun

  print,'Total time consuming is:',systime(1)-start_time,'seconds'
  print,'This program is successful performed!'

end