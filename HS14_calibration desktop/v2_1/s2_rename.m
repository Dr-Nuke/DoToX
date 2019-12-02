% renames the files, copies them to \processd & fills in the xls log file
     
if any(itr_cas==1)
    s2_rename_emp
end

if any(itr_cas==2)
    s2_rename_d2o
end

if any(itr_cas==3)
    s2_rename_cl3
end

if any(itr_cas==4)
    s2_rename_dc_
end

if any(itr_cas==5)
    s2_rename_ob_
end

xls={emp_xls,d2o_xls,cl3_xls};
save('xls_file','xls');


; %1

; %2

; %3

; %4

; %5

; %6

; %7

; %8

; %9

; %10