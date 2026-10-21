PRO Detrend_001

  start_time=systime(1)

  infile_Tmp ='E:\Work22\09_Regress2\D2\RH_VEG_LR.dat'
;  infile_Pre ='D:\MOD13_Yellow_River\MCD12Q1\output2\Pre_2001_2020_02.dat'
;  infile_NDVI ='D:\MOD13_Yellow_River\MOD13A3\Detrend\NDVI_2001_2020.dat'
  infile_Tmp_K='E:\Work22\09_Regress2\D2\Trend\RH_VEG_LR_Trend.dat'
;  infile_Pre_K='D:\MOD13_Yellow_River\MOD13A3\Detrend\Pre.dat'
;  infile_NDVI_K='D:\MOD13_Yellow_River\MOD13A3\Detrend\NDVI.dat'
  
  column = 320
  row = 160
  band1 = 13
  band2 = 1

  Tmp_data = read_binary(infile_Tmp,data_dims = [column,row,band1],data_type = 4)
;  Pre_data = read_binary(infile_Pre,data_dims = [column,row,band1],data_type = 4)
;  NDVI_data = read_binary(infile_NDVI,data_dims = [column,row,band1],data_type = 4)
  Tmp_K_data = read_binary(infile_Tmp_K,data_dims = [column,row,band2],data_type = 4)
;  Pre_K_data = read_binary(infile_Pre_K,data_dims = [column,row,band2],data_type = 4)
;  NDVI_K_data = read_binary(infile_NDVI_K,data_dims = [column,row,band2],data_type = 4)
  ;result_SM = make_array(column,row,type = 4)

  ;创建输出变量
  Detrend1=make_array(column,row,band1, type=4, value=0.0);动态创建数组
  Detrend2=make_array(2795,1152,20, type=4, value=0.0);动态创建数组
  Detrend3=make_array(2795,1152,20, type=4, value=0.0);动态创建数组
  
    for i = 0,column-1 do begin
      for j = 0,row-1 do begin
        for k = 0,13-1 do begin
        x1 = reform(Tmp_data(i,j,(0+k):(0+k)))
        x2 = reform(Tmp_K_data(i,j,0))
;        y1 = reform(Pre_data(i,j,(0+k):(0+k)))
;        y2 = reform(Pre_K_data(i,j,0))
;        z1 = reform(NDVI_data(i,j,(0+k):(0+k)))
;        z2 = reform(NDVI_K_data(i,j,0))
        
        Detrend1[i,j,k]=x1-x2*(k+1)
;        Detrend2[i,j,k]=y1-y2*(k+1)
;        Detrend3[i,j,k]=z1-z2*(k+1)
        
      endfor
    endfor
  endfor
  
  out_file1=strcompress('E:\Work22\08_Regress\Detrend02\'+'RH_VEG_LR_Detrend'+'.dat',/remove_all)
  openw,out_lun,out_file1,/get_lun
  writeu,out_lun,Detrend1
  free_lun,out_lun
  
;  out_file2=strcompress('D:\MOD13_Yellow_River\MOD13A3\Detrend\'+'Pre_Detrend'+'.dat',/remove_all)
;  openw,out_lun,out_file2,/get_lun
;  writeu,out_lun,Detrend2
;  free_lun,out_lun
;  
;  out_file3=strcompress('D:\MOD13_Yellow_River\MOD13A3\Detrend\'+'NDVI_Detrend'+'.dat',/remove_all)
;  openw,out_lun,out_file3,/get_lun
;  writeu,out_lun,Detrend3
;  free_lun,out_lun

  print,'Total time consuming is:',systime(1)-start_time,'seconds!'

  print,'finished!'
 END