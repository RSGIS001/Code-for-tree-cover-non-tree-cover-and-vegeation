PRO R_correlation

  start_time=systime(1)

  ;infile_WUE = 'F:\SM_VPD_0202\20210202\Global_Month\tmp_1982_2015_Monthly.dat'
  infile_NDVI = 'E:\Work3\04_Basic_Data3\Trend\EI2.dat'
  ;infile_LAI = 'F:\SM_VPD_0202\20210202\Global_Month\VPD_1982_2015_Monthly.dat'
  infile_Prec = 'E:\Work3\04_Basic_Data3\Trend\ET2.dat'
 ; infile_ET = 'H:\20210717\CRU_Pre\S6\VPD_AI_VEG_Anomaly.dat'

  column = 1440
  row = 720
  band = 11

  ;WUE_data = read_binary(infile_WUE,data_dims = [column,row,band],data_type = 4)
  NDVI_data = read_binary(infile_NDVI,data_dims = [column,row,band],data_type = 4)
;  LAI_data = read_binary(infile_LAI,data_dims = [column,row,band],data_type = 4)
  Prec_data = read_binary(infile_Prec,data_dims = [column,row,band],data_type = 4)
;  ET_data = read_binary(infile_ET,data_dims = [column,row,band],data_type = 4)

  ;定义结果数组
  result_1 = make_array(column,row,2, type = 4,value = 0.0)
  result_2 = make_array(column,row,2, type = 4,value = 0.0)
  result_3 = make_array(column,row,2, type = 4,value = 0.0)
  ;result_LAI = make_array(column,row,type = 5)
  ;result_VPD = make_array(column,row,type = 4)
  ;result_TMP = make_array(column,row,type = 4)
  ;result_max1 = make_array(column,row,type = 4,value = 9999.0)
  ;result_max2 = make_array(column,row,type = 4,value = 9999.0)
  
  ;calculate partial correlation
  for i = 0,column-1 do begin
    for j = 0,row-1 do begin
      x = reform(NDVI_data(i,j,*))
      y1 = reform(Prec_data(i,j,*))
  ;    y2 = reform(ET_data(i,j,*))
      
      y11=(y1-mean(y1))/stddev(y1)
      
          result_1[i,j,*] = R_correlate(x,y1) 
 ;         result_2[i,j,*] = R_correlate(x,y2)
 ;         result_3[i,j,*] = R_correlate(x,y2) 

    endfor
  endfor

  ;写出结果
  outfile1 = strcompress('E:\Work3\04_Basic_Data3\Trend_Sig\'+'EI_ET_025_CORREL'+'.dat',/remove_all);
  openw,out_lun,outfile1,/get_lun
  writeu,out_lun,result_1
  free_lun,out_lun
;  outfile2 = strcompress('H:\20210717\CRU_Pre\S6\CORREL\'+'NDVI_VPD_Sig'+'.dat',/remove_all);
;  openw,out_lun,outfile2,/get_lun
;  writeu,out_lun,result_2
;  free_lun,out_lun
;  outfile3 = strcompress('F:\SM_VPD_0202\20210202\Output02\05_Trend_Mapping\Point2022\'+'SM_T'+'.dat',/remove_all);
;  openw,out_lun,outfile3,/get_lun
;  writeu,out_lun,result_3
;  free_lun,out_lun
;  outfile = strcompress('F:\SM_VPD_China\Model\Correlation\'+'NDVI_VPD_1982_2015'+'.dat',/remove_all);
;  openw,out_lun,outfile,/get_lun
;  writeu,out_lun,result_LAI
;  free_lun,out_lun
;  outfile = strcompress('E:\yu\重采样结果\'+'NDVI_TmpMAX_2000_2015'+'.dat',/remove_all);
;  openw,out_lun,outfile,/get_lun
;  writeu,out_lun,result_VPD
;  free_lun,out_lun
;  outfile = strcompress('E:\yu\重采样结果\'+'NDVI_TmpMIN_2000_2015'+'.dat',/remove_all);
;  openw,out_lun,outfile,/get_lun
;  writeu,out_lun,result_TMP
;  free_lun,out_lun
;  outfile = strcompress('F:\10-partial_correlation\'+'MAX_CC_02_11'+'.dat',/remove_all);
;  openw,out_lun,outfile,/get_lun
;  writeu,out_lun,result_max1
;  free_lun,out_lun
;  outfile = strcompress('F:\10-partial_correlation\'+'MAX2_92_01'+'.dat',/remove_all);
;  openw,out_lun,outfile,/get_lun
;  writeu,out_lun,result_max2
;  free_lun,out_lun  

  print,'Total time consuming is:',systime(1)-start_time,'seconds!'
  print,'finished!'

end