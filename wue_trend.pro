;+
;This program is used to calculate mean and trend of wue
;
;method: Sen slope
;
;Date: 2017-3-27
;+
FUNCTION mkstrend,x,undef=undef
  
  if total(finite(x,/nan)) eq 0 then begin
    
    x = float(x)
    
    nx = n_elements(x)
    
    nx1 = nx - 1.
    
    n = nx*(nx - 1)/2.  ; The number of elements in d
    
    d = fltarr(n)
    
    m = 0.
    
    for i = 0,nx-1 do begin
      
      for j = i+1,nx-1 do begin
        
        d(m) = x(j) - x(i)
        
        m = m + 1
        
      endfor
      
    endfor
    
    for i = 0L,n-1 do begin
      
      if d(i) lt 0 then d(i) = -1.
      
      if d(i) eq 0 then d(i) = 0.
      
      if d(i) gt 0 then d(i) = 1.
      
    endfor
    
    s = total(d)
    
    U = x(uniq(x(sort(x))))
    
    Corr = 0 ;Correction for tied observations(equal value)
    
    for y = 0,n_elements(U) - 1 do begin
      
      find = where(x eq U(y))
      
      uj = n_elements(find)
      
      Corr = Corr + uj * (uj - 1) * (2 * uj + 5)
      
    endfor
    
    Vs = (nx * (nx - 1.) * (2 * nx + 5.) - Corr)/18.
    
    if s gt 0. then z = (s - 1)/sqrt(Vs)
    
    if s lt 0. then z = (s + 1)/sqrt(Vs)
    
    if s eq 0. then z = 0.
    
    nor = gauss_cvf(0.025) ; Prob at 95%(two-side)  [0,1]
    
    Sn = fltarr(n)
    
    m = 0.
    
    for i = 0,nx1-1 do begin
      
      for j = i + 1,nx-1 do begin
        
        Sn(m) = (x(i)-x(j))/(i-j)
        
        m = m + 1
        
      endfor
      
    endfor
    
    Snsorted = Sn(sort(Sn))
    
    m = float(fix(n/2.))
    
    if abs(z) lt nor then begin
      
      slope = undef
      
    endif else begin
      
      if 2 * m eq n then slope = 0.5 * (Snsorted(m) + Snsorted(m + 1))
      
      if 2 * m + 1. eq n then slope = Snsorted(m + 1)
      
    endelse  
    
  endif else begin
    
    slope = undef
    
  endelse
  
  return,slope
  
END

;+
;
;The following is the main program
;
;+

;pro wue_trend

  start_time=systime(1)
  
  infile_wue ='E:\Work3\04_Basic_Data3\Trend\EI2.dat'
  
 column = 1440
 
 row = 720
 
 band = 11
  
  WUE_data = make_array(column,row,band,type = 4)
  ;WUE_data = read_binary(column,row,band,type = 4,value = 32767.0)

  ;璇�UE�版�

  openr,in_lun,infile_wue,/get_lun

  readu,in_lun,WUE_data

  free_lun,in_lun
  
  ;Define result array
  
  ;result_mean = make_array(column,row,type = 4,value = 9999.0)  ; ���
  
  result_mk = make_array(column,row,type = 4)  ; Sen瓒��
  
  ;result_slope = make_array(column,row,type = 4,value = 9999.0)  ; 绾挎�瓒��
  
  ;result_Spr = make_array(column,row,type = 4,value = 9999.0)   ; Spearman绉╃��崇郴��
  
  ;result_F = make_array(column,row,type = 4,value = 9999.0)   ; F妫��
  
  ;result_T = make_array(column,row,type = 4,value = 9999.0)   ; T妫��
  
  ;���������绠�
  
  for ii = 0,column-1 do begin
    
    for jj = 0,row-1 do begin
      
      x = reform(WUE_data(ii,jj,*))
      
      index = where(x lt 0)
      
      if index[0] eq -1 then begin
      
        ;mean wue
      
        ;result_mean[ii,jj] = mean(x)   

        ;MK trend

        result_mk[ii,jj] = mkstrend(x,undef = 0)
      
        ;Linear trend
      
        ;time = indgen(band)
      
        ;coeff = LINFIT(time,x,YFIT = x_predicted)
      
        ;result_slope[ii,jj] = coeff[1]
      
        ;F test
      
        ;if total((x-x_predicted)^2) ne 0.0 then begin
        
          ;result_F[ii,jj] = total((x_predicted - mean(x))^2)/(total((x-x_predicted)^2)/(band -2))
        
        ;endif else begin
        
          ;result_F[ii,jj] = 0.0
        
        ;endelse
;      
;        ;Spearman test
;      
;        y = where(x eq x[0],count)
;      
;        if count ne band then begin
;        
;          y = R_CORRELATE(x,time)
;        
;          result_Spr[ii,jj] = y[0]
;        
;        endif else begin
;        
;          result_Spr[ii,jj] = 0.0
;        
;        endelse
;      
;        ;t test
;      
;        T[ii,jj] = abs(result_Spr[ii,jj])*sqrt((band - 2)/(1 - (result_Spr[ii,jj])^2))  
      
      endif else begin

        ;result_mean[ii,jj] = float(9999)
        
        result_mk[ii,jj] = float(0.0)
        
        ;result_slope[ii,jj] = float(9999)
        
        ;result_F[ii,jj] = float(9999)
;        
;        result_Spr[ii,jj] = float(9999)
;        
;        result_t[ii,jj] = float(9999)

      endelse 
      
    endfor
       
  endfor
  
  ;���缁��
  
  ;Mean WUE

;  outfile_mean = strcompress('G:\WUE_00_14\'+'result_mean'+'.dat',/remove_all);
;
;  openw,out_lun,outfile_mean,/get_lun
;
;  writeu,out_lun,result_mean
;
;  free_lun,out_lun 
   
  ;MK trend
  
  outfile_mk = strcompress('E:\Work3\04_Basic_Data3\Trend_Sig\'+'EI2_Trend_005'+'.dat',/remove_all);
  
  openw,out_lun,outfile_mk,/get_lun

  writeu,out_lun,result_mk

  free_lun,out_lun

  ;linear trend

;  outfile_linear = strcompress('G:\WUE_00_14\'+'result_slope'+'.dat',/remove_all);
;
;  openw,out_lun,outfile_linear,/get_lun
;
;  writeu,out_lun,result_slope
;
;  free_lun,out_lun
  
  ;F test

;  outfile_f = strcompress('G:\WUE_00_14\'+'result_F'+'.dat',/remove_all);
;
;  openw,out_lun,outfile_f,/get_lun
;
;  writeu,out_lun,result_F
;
;  free_lun,out_lun
;  
;  ;Spearman test
;
;  outfile_Spr = strcompress('G:\WUE_00_14\'+'result_Spr'+'.dat',/remove_all);
;
;  openw,out_lun,outfile_Spr,/get_lun
;
;  writeu,out_lun,result_Spr
;
;  free_lun,out_lun
;  
;  ;t test
;
;  outfile_t = strcompress('G:\WUE_00_14\'+'result_T'+'.dat',/remove_all);
;
;  openw,out_lun,outfile_t,/get_lun
;
;  writeu,out_lun,result_T
;
;  free_lun,out_lun
  
  
  print,'Total time consuming is:',systime(1)-start_time,'seconds!'

  print,'finished!'
  
end