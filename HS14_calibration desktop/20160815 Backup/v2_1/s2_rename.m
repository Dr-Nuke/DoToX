% renames the files, copies them to \processd & fills in the xls log file

ind_copy=1;
        if ind_copy==0
            disp('not copying files. see s2_rename.m')
        end

s2_rename_emp
s2_rename_d2o
s2_rename_cl3
s2_rename_DC

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