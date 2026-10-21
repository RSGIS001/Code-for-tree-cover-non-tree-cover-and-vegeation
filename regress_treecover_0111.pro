PRO Regress_TreeCover_0111


  start_time=systime(1)
  

  
  
  infile_TC =  'E:\Work22\08_Regress\Detrend02\TC_BG_EC-Earth3-CC_alog_Detrend.dat'
  infile_NTC = 'E:\Work22\08_Regress\Detrend02\NTC_BG_EC-Earth3-CC_alog_Detrend.dat'
  infile_VEG =  'E:\Work22\08_Regress\Detrend02\RH_CC_Detrend.dat'
  infile_AI = 'E:\Work22\08_Regress\Detrend02\AI_EC-Earth3-CC_Detrend.dat'
  infile_Ta = 'E:\Work22\08_Regress\Detrend02\TAS_EC-Earth3-CC_Detrend.dat'
  infile_SR = 'E:\Work22\08_Regress\Detrend02\RSS_EC-Earth3-CC_Detrend.dat'
  infile_WS = 'E:\Work22\08_Regress\Detrend02\WS_EC-Earth3-CC_Detrend.dat'
  infile_CO2 = 'E:\Work22\08_Regress\Detrend02\CO2_Detrend_512_256.dat'
  
  
  column = 512
  row = 256
  band = 13

  TC_data = read_binary(infile_TC,data_dims = [column,row,band],data_type = 4)
  NTC_data = read_binary(infile_NTC,data_dims = [column,row,band],data_type = 4)
  VEG_data = read_binary(infile_VEG,data_dims = [column,row,band],data_type = 4)
  AI_data = read_binary(infile_AI,data_dims = [column,row,band],data_type = 4)
  Ta_data = read_binary(infile_Ta,data_dims = [column,row,band],data_type = 4)
  SR_data = read_binary(infile_SR,data_dims = [column,row,band],data_type = 4)
  WS_data = read_binary(infile_WS,data_dims = [column,row,band],data_type = 4)
  CO2_data = read_binary(infile_CO2,data_dims = [column,row,band],data_type = 4)
  
  ;定义结果数组
  result_1 = make_array(column,row,7, type = 4,value = 0.0)
  result_2 = make_array(column,row,7, type = 4,value = 0.0)
  result_3 = make_array(column,row,2, type = 4,value = 0.0)
  result_4 = make_array(column,row,13, type = 4,value = 0.0)
  result_5 = make_array(column,row,13, type = 4,value = 0.0)
  result_6 = make_array(column,row,7, type = 4,value = 0.0)
  
  
  for i = 0,column-1 do begin
    for j = 0,row-1 do begin
  
      Y10 = reform(VEG_data(i,j,*))
      X10 = reform(TC_data(i,j,*))
      X20 = reform(NTC_data(i,j,*))
      X30 = reform(AI_data(i,j,*))
      X40 = reform(Ta_data(i,j,*))
      X50 = reform(SR_data(i,j,*))
      X60 = reform(WS_data(i,j,*))
      X70 = reform(CO2_data(i,j,*))
;      X80 = reform(SM_data(i,j,*))
      
      Y = (Y10-mean(Y10))/stddev(Y10)
      X1 = (X10-mean(X10))/stddev(X10)
      X2 = (X20-mean(X20))/stddev(X20)
      X3 = (X30-mean(X30))/stddev(X30)
      X4 = (X40-mean(X40))/stddev(X40)
      X5 = (X50-mean(X50))/stddev(X50)
      X6 = (X60-mean(X60))/stddev(X60)
      X7 = (X70-mean(X70))/stddev(X70)
;      X8 = (X80-mean(X80))/stddev(X80)
      
      Y100 = stddev(Y10)
      X100 = stddev(X10)
      X200 = stddev(X20)
      X300 = stddev(X30)
      X400 = stddev(X40)
      X500 = stddev(X50)
      X600 = stddev(X60)
      X700 = stddev(X70)
 ;     X800 = stddev(X80)
      
      
      
      ;XX=[TRANSPOSE(X1),TRANSPOSE(X2),TRANSPOSE(X3)]
      XX0=[TRANSPOSE(X1),TRANSPOSE(X2),TRANSPOSE(X3),TRANSPOSE(X4),TRANSPOSE(X5),TRANSPOSE(X6),TRANSPOSE(X7)]
      XX1=[TRANSPOSE(X2),TRANSPOSE(X3),TRANSPOSE(X4),TRANSPOSE(X5),TRANSPOSE(X6),TRANSPOSE(X7)]
      XX2=[TRANSPOSE(X1),TRANSPOSE(X3),TRANSPOSE(X4),TRANSPOSE(X5),TRANSPOSE(X6),TRANSPOSE(X7)]
      XX3=[TRANSPOSE(X1),TRANSPOSE(X2),TRANSPOSE(X4),TRANSPOSE(X5),TRANSPOSE(X6),TRANSPOSE(X7)]
      XX4=[TRANSPOSE(X1),TRANSPOSE(X2),TRANSPOSE(X3),TRANSPOSE(X5),TRANSPOSE(X6),TRANSPOSE(X7)]
      XX5=[TRANSPOSE(X1),TRANSPOSE(X2),TRANSPOSE(X3),TRANSPOSE(X4),TRANSPOSE(X6),TRANSPOSE(X7)]
      XX6=[TRANSPOSE(X1),TRANSPOSE(X2),TRANSPOSE(X3),TRANSPOSE(X4),TRANSPOSE(X5),TRANSPOSE(X7)]
      XX7=[TRANSPOSE(X1),TRANSPOSE(X2),TRANSPOSE(X3),TRANSPOSE(X4),TRANSPOSE(X5),TRANSPOSE(X6)]
;      XX8=[TRANSPOSE(X1),TRANSPOSE(X2),TRANSPOSE(X3),TRANSPOSE(X4),TRANSPOSE(X5),TRANSPOSE(X6)]
      
;      index = where(Y eq 0)
;      if index[0] eq -1 then begin
;      endif
      

      
      ;初始化高斯误差
      measure_errors = REPLICATE(0.5, N_ELEMENTS(Y))
      ;多元回归
;      result0 = REGRESS(XX, X7, SIGMA=sigma0, CONST=const0, MEASURE_ERRORS=measure_errors0, $
;        STATUS=status0, FTEST=ftest0, CORRELATION=correlation0, MCORRELATION=mcorrelation0)      
;      P=mcorrelation0^2
         
      result00 = REGRESS(XX0, Y, SIGMA=sigma00, CONST=const00, MEASURE_ERRORS=measure_errors00, $
        STATUS=status00, FTEST=ftest00, CORRELATION=correlation00, MCORRELATION=mcorrelation00)
      result11 = REGRESS(XX1, Y, SIGMA=sigma11, CONST=const11, MEASURE_ERRORS=measure_errors11, $
        STATUS=status11, FTEST=ftest11, CORRELATION=correlation11, MCORRELATION=mcorrelation11)
      result22 = REGRESS(XX2, Y, SIGMA=sigma22, CONST=const22, MEASURE_ERRORS=measure_errors22, $
        STATUS=status22, FTEST=ftest22, CORRELATION=correlation22, MCORRELATION=mcorrelation22)
      result33 = REGRESS(XX3, Y, SIGMA=sigma33, CONST=const33, MEASURE_ERRORS=measure_errors33, $
        STATUS=status33, FTEST=ftest33, CORRELATION=correlation33, MCORRELATION=mcorrelation33)
      result44 = REGRESS(XX4, Y, SIGMA=sigma44, CONST=const44, MEASURE_ERRORS=measure_errors44, $
        STATUS=status44, FTEST=ftest44, CORRELATION=correlation44, MCORRELATION=mcorrelation44)
      result55 = REGRESS(XX5, Y, SIGMA=sigma55, CONST=const55, MEASURE_ERRORS=measure_errors55, $
        STATUS=status55, FTEST=ftest55, CORRELATION=correlation55, MCORRELATION=mcorrelation55)
      result66 = REGRESS(XX6, Y, SIGMA=sigma66, CONST=const66, MEASURE_ERRORS=measure_errors66, $
        STATUS=status66, FTEST=ftest66, CORRELATION=correlation66, MCORRELATION=mcorrelation66)
      result77 = REGRESS(XX7, Y, SIGMA=sigma77, CONST=const77, MEASURE_ERRORS=measure_errors77, $
        STATUS=status77, FTEST=ftest77, CORRELATION=correlation77, MCORRELATION=mcorrelation77)
;      result88 = REGRESS(XX8, Y, SIGMA=sigma88, CONST=const88, MEASURE_ERRORS=measure_errors88, $
;        STATUS=status88, FTEST=ftest88, CORRELATION=correlation88, MCORRELATION=mcorrelation88)
      
      result111 = REGRESS(XX1, X1, MCORRELATION=mcorrelation111)
      result222 = REGRESS(XX2, X2, MCORRELATION=mcorrelation222)
      result333 = REGRESS(XX3, X3, MCORRELATION=mcorrelation333)
      result444 = REGRESS(XX4, X4, MCORRELATION=mcorrelation444)
      result555 = REGRESS(XX5, X5, MCORRELATION=mcorrelation555)
      result666 = REGRESS(XX6, X6, MCORRELATION=mcorrelation666)
      result777 = REGRESS(XX7, X7, MCORRELATION=mcorrelation777)
      ;result888 = REGRESS(XX7, X8, MCORRELATION=mcorrelation888)
      
      
      
      d0=const00
      e0=correlation00
      TC_Partial_R2= (result00[0,0]/sigma00[0])^2/((result00[0,0]/sigma00[0])^2+5)
      NTC_Partial_R2=(result00[0,1]/sigma00[1])^2/((result00[0,1]/sigma00[1])^2+5)
      X3_Partial_R2=(result00[0,2]/sigma00[2])^2/((result00[0,2]/sigma00[2])^2+5)
      X4_Partial_R2=(result00[0,3]/sigma00[3])^2/((result00[0,3]/sigma00[3])^2+5)
      X5_Partial_R2=(result00[0,4]/sigma00[4])^2/((result00[0,4]/sigma00[4])^2+5)
      X6_Partial_R2=(result00[0,5]/sigma00[5])^2/((result00[0,5]/sigma00[5])^2+5)
      X7_Partial_R2=(result00[0,6]/sigma00[6])^2/((result00[0,6]/sigma00[6])^2+5)
     ; X8_Partial_R2=(result00[0,7]/sigma00[7])^2/((result00[0,7]/sigma00[7])^2+5)
      
      
      
      R2=mcorrelation00^2
      pvalue = 1-F_PDF(ftest00, 7, 5)
      
      Result_1[i,j,0]= Y100*result00[0,0]/X100
      Result_1[i,j,1]= Y100*result00[0,1]/X200
      Result_1[i,j,2]= Y100*result00[0,2]/X300
      Result_1[i,j,3]= Y100*result00[0,3]/X400
      Result_1[i,j,4]= Y100*result00[0,4]/X500
      Result_1[i,j,5]= Y100*result00[0,5]/X600
      Result_1[i,j,6]= Y100*result00[0,6]/X700
      ;Result_1[i,j,7]= Y100*result00[0,7]/X800
      
      Result_2[i,j,0]= TC_Partial_R2
      Result_2[i,j,1]= NTC_Partial_R2
      Result_2[i,j,2]= X3_Partial_R2
      Result_2[i,j,3]= X4_Partial_R2
      Result_2[i,j,4]= X5_Partial_R2
      Result_2[i,j,5]= X6_Partial_R2
      Result_2[i,j,6]= X7_Partial_R2
      ;Result_2[i,j,7]= X8_Partial_R2
      
      Result_3[i,j,0]= R2
      Result_3[i,j,1]= pvalue
      Result_4[i,j,*]= d0+result00[0,0]*X1+result00[0,1]*X2+result00[0,2]*X3+ $
        result00[0,3]*X4+result00[0,4]*X5+result00[0,5]*X6+result00[0,6]*X7
      Result_5[i,j,*]= Y-(d0+result00[0,0]*X1+result00[0,1]*X2+result00[0,2]*X3+ $
        result00[0,3]*X4+result00[0,4]*X5+result00[0,5]*X6+result00[0,6]*X7)
      
      Result_6[i,j,0]= 1/(1-mcorrelation111^2)
      Result_6[i,j,1]= 1/(1-mcorrelation222^2)
      Result_6[i,j,2]= 1/(1-mcorrelation333^2)
      Result_6[i,j,3]= 1/(1-mcorrelation444^2)
      Result_6[i,j,4]= 1/(1-mcorrelation555^2)
      Result_6[i,j,5]= 1/(1-mcorrelation666^2)
      Result_6[i,j,6]= 1/(1-mcorrelation777^2)
      ;Result_6[i,j,7]= 1/(1-mcorrelation888^2)
  
  
  
      
    endfor
  endfor
  
  ;写出结果
  outfile1 = strcompress('E:\Work22\08_Regress\16EC-Earth3-CC\RH\'+'TC_NTC_TCNTC_Sensitivity'+'.dat',/remove_all);
  openw,out_lun,outfile1,/get_lun
  writeu,out_lun,result_1
  free_lun,out_lun
  
  outfile2 = strcompress('E:\Work22\08_Regress\16EC-Earth3-CC\RH\'+'TC_NTC_Partial_R2'+'.dat',/remove_all);
  openw,out_lun,outfile2,/get_lun
  writeu,out_lun,result_2
  free_lun,out_lun
  
  outfile3 = strcompress('E:\Work22\08_Regress\16EC-Earth3-CC\RH\'+'TC_NTC_R2_Pvalue'+'.dat',/remove_all);
  openw,out_lun,outfile3,/get_lun
  writeu,out_lun,result_3
  free_lun,out_lun
  
  outfile4 = strcompress('E:\Work22\08_Regress\16EC-Earth3-CC\RH\'+'VEG_Sim'+'.dat',/remove_all);
  openw,out_lun,outfile4,/get_lun
  writeu,out_lun,result_4
  free_lun,out_lun
  
  outfile5 = strcompress('E:\Work22\08_Regress\16EC-Earth3-CC\RH\'+'VEG_Sim_Residual'+'.dat',/remove_all);
  openw,out_lun,outfile5,/get_lun
  writeu,out_lun,result_5
  free_lun,out_lun
  
  outfile6 = strcompress('E:\Work22\08_Regress\16EC-Earth3-CC\RH\'+'VIF'+'.dat',/remove_all);
  openw,out_lun,outfile6,/get_lun
  writeu,out_lun,result_6
  free_lun,out_lun
  
  
  print,'Total time consuming is:',systime(1)-start_time,'seconds!'
  print,'finished!'


END