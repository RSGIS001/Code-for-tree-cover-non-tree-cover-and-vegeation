PRO Regress_TreeCover_02


  start_time=systime(1)

  infile_TC =  'E:\Work22\02_VCF\050\TC_BG_050_alog.dat'
  infile_NTC = 'E:\Work22\02_VCF\050\NTC_BG_050_alog.dat'
  infile_TCNTC = 'E:\Work22\02_VCF\050\alog\TC_NTC_050_alog_Centered.dat'
  infile_VEG =  'E:\Work22\04_LAI\LAI3g_VEG.dat'
  infile_AI = 'E:\Work22\03_CRU\CRU\AI_2001_2013.dat'
  infile_PET = 'E:\Work22\03_CRU\CRU\PET_2001_2013.dat'
  infile_Ta = 'E:\Work22\03_CRU\CRU\Ta_2001_2013.dat'
  infile_Cloud = 'E:\Work22\03_CRU\CRU\Cloud_2001_2013.dat'
  infile_SM = 'E:\Work22\03_CRU\ERA\SM_2001_2013.dat'
  infile_VPD = 'E:\Work22\03_CRU\CRU\VPD_2001_2013.dat'
  infile_CO2 = 'E:\Work22\03_CRU\CRU\CO2_2001_2013.dat'
  infile_SR = 'E:\Work22\03_CRU\ERA\SR_2001_2013.dat'
  infile_WS = 'E:\Work22\03_CRU\ERA\WS_2001_2013.dat'

  column = 720
  row = 360
  band = 13

  TC_data = read_binary(infile_TC,data_dims = [column,row,band],data_type = 4)
  NTC_data = read_binary(infile_NTC,data_dims = [column,row,band],data_type = 4)
  VEG_data = read_binary(infile_VEG,data_dims = [column,row,band],data_type = 4)
  AI_data = read_binary(infile_AI,data_dims = [column,row,band],data_type = 4)
  Ta_data = read_binary(infile_Ta,data_dims = [column,row,band],data_type = 4)
  Cloud_data = read_binary(infile_Cloud,data_dims = [column,row,band],data_type = 4)
  SR_data = read_binary(infile_SR,data_dims = [column,row,band],data_type = 4)
  CO2_data = read_binary(infile_CO2,data_dims = [column,row,band],data_type = 4)
  PET_data = read_binary(infile_PET,data_dims = [column,row,band],data_type = 4)
  WS_data = read_binary(infile_WS,data_dims = [column,row,band],data_type = 4)

  ;定义结果数组
  result_1 = make_array(column,row,7, type = 4,value = 0.0)
  result_2 = make_array(column,row,2, type = 4,value = 0.0)
  result_3 = make_array(column,row,7, type = 4,value = 0.0)
  result_4 = make_array(column,row,7, type = 4,value = 0.0)
  result_5 = make_array(column,row,13, type = 4,value = 0.0)
  result_6 = make_array(column,row,13, type = 4,value = 0.0)


  for i = 0,column-1 do begin
    for j = 0,row-1 do begin

      Y = reform(VEG_data(i,j,*))
      X1 = reform(TC_data(i,j,*))
      X2 = reform(Cloud_data(i,j,*))
      X3 = reform(PET_data(i,j,*))
      X4 = reform(Ta_data(i,j,*))
      X5 = reform(SR_data(i,j,*))
      X6 = reform(WS_data(i,j,*))
      X7 = reform(CO2_data(i,j,*))
      X8 = reform(Cloud_data(i,j,*))
      
      
      
      index = where(Y eq 0)
      if index[0] eq -1 then begin
      endif

      XX0=[TRANSPOSE(X1),TRANSPOSE(X2),TRANSPOSE(X3),TRANSPOSE(X4),TRANSPOSE(X5),TRANSPOSE(X6),TRANSPOSE(X7)]
     
      
      PRINCOMP2,in_data = XX0,loadings = loadings,scores= scores,latent= latent

      a0=latent
      b0=scores
      c0=loadings

      result_1[i,j,*] = a0
      
      XX1= b0[0,*]
      XX2= b0[1,*]
      XX3= b0[2,*]
      XX= [XX1,XX2,XX3]

      ;初始化高斯误差
      measure_errors = REPLICATE(0.5, N_ELEMENTS(Y))

      ;多元回归

      result00 = REGRESS(XX, Y, SIGMA=sigma00, CONST=const00, MEASURE_ERRORS=measure_errors00, $
        STATUS=status00, FTEST=ftest00, CORRELATION=correlation00, MCORRELATION=mcorrelation00)
      



      d0=const00
      e0=correlation00
      
      R2=mcorrelation00^2
      pvalue = 1-F_PDF(ftest00, 3, 9)

      Result_2[i,j,0]= R2
      Result_2[i,j,1]= pvalue     
      Result_3[i,j,*]=result00[0,0]*loadings[0,*]+result00[0,1]*loadings[1,*]+result00[0,2]*loadings[2,*]
      Result_4[i,j,*]=abs(result00[0,0]*loadings[0,*])+abs(result00[0,1]*loadings[1,*])+abs(result00[0,2]*loadings[2,*])     
      Result_5[i,j,*]=d0+result00[0,0]*XX1+result00[0,1]*XX2+result00[0,2]*XX3
      Result_6[i,j,*]= Y-(d0+result00[0,0]*XX1+result00[0,1]*XX2+result00[0,2]*XX3)

      




    endfor
  endfor

  ;写出结果
  outfile1 = strcompress('E:\Work22\08_Regress\1LAI4g\PCA\'+'PCs_STDDEV'+'.dat',/remove_all);
  openw,out_lun,outfile1,/get_lun
  writeu,out_lun,result_1
  free_lun,out_lun

  outfile2 = strcompress('E:\Work22\08_Regress\1LAI4g\PCA\'+'TC_NTC_R2_Pvalue'+'.dat',/remove_all);
  openw,out_lun,outfile2,/get_lun
  writeu,out_lun,result_2
  free_lun,out_lun

  outfile3 = strcompress('E:\Work22\08_Regress\1LAI4g\PCA\'+'VEG_Sensitivity'+'.dat',/remove_all);
  openw,out_lun,outfile3,/get_lun
  writeu,out_lun,result_3
  free_lun,out_lun

  outfile4 = strcompress('E:\Work22\08_Regress\1LAI4g\PCA\'+'VEG_Dominant_Factor'+'.dat',/remove_all);
  openw,out_lun,outfile4,/get_lun
  writeu,out_lun,result_4
  free_lun,out_lun

  outfile5 = strcompress('E:\Work22\08_Regress\1LAI4g\PCA\'+'VEG_Sim'+'.dat',/remove_all);
  openw,out_lun,outfile5,/get_lun
  writeu,out_lun,result_5
  free_lun,out_lun

  outfile6 = strcompress('E:\Work22\08_Regress\2LAI3g\'+'VEG_Sim_Residual'+'.dat',/remove_all);
  openw,out_lun,outfile6,/get_lun
  writeu,out_lun,result_6
  free_lun,out_lun


  print,'Total time consuming is:',systime(1)-start_time,'seconds!'
  print,'finished!'


END