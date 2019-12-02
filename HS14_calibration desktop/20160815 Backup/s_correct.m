
    % DC correction
val_smooth=0.1;    
    
if ind_debug==1, %see beginning total.m

    d2o_corr     =f_DC(d2o,DC);
    d2o_180_corr =f_DC(d2o_180,DC);
    d2o_OB_corr  =f_DC(d2o_OB,DC);

    cl3_corr     =f_DC(cl3,DC);
    cl3_180_corr =f_DC(cl3_180,DC);
    cl3_OB_corr  =f_DC(cl3_OB,DC);

    empty_corr     =f_DC(empty,DC);
    empty_180_corr =f_DC(empty_180,DC);
    empty_OB_corr  =f_DC(empty_OB,DC);


    % flat field correction
    d2o_corr=f_ffcorrect(d2o_OB_corr,d2o_corr);
    d2o_180_corr=f_ffcorrect(d2o_OB_corr,d2o_180_corr);

    cl3_corr=f_ffcorrect(d2o_OB_corr,cl3_corr);
    cl3_180_corr=f_ffcorrect(cl3_OB_corr,cl3_180_corr);

    empty_corr=f_ffcorrect(empty_OB_corr,empty_corr);
    empty_180_corr=f_ffcorrect(empty_OB_corr,empty_180_corr);



    %dose correction
    d2o_corr=f_doseCorrect(d2o_corr,xdose,ydose);
    d2o_180_corr=f_doseCorrect(d2o_180_corr,xdose,ydose);
    d2o_OB_corr=f_doseCorrect(d2o_OB,xdose,ydose);

    cl3_corr=f_doseCorrect(cl3_corr,xdose,ydose);
    cl3_180_corr=f_doseCorrect(cl3_180_corr,xdose,ydose);
    cl3_OB_corr=f_doseCorrect(cl3_OB,xdose,ydose);

    empty_corr=f_doseCorrect(empty_corr,xdose,ydose);
    empty_180_corr=f_doseCorrect(empty_180_corr,xdose,ydose);
    empty_OB_corr=f_doseCorrect(empty_OB,xdose,ydose);
    
    
    %% experimental: ring artefact correction
    d2o_corr=f_ringCorr(d2o_corr,val_smooth);
    d2o_180_corr=f_ringCorr(d2o_180_corr,val_smooth);
    d2o_OB_corr=f_ringCorr(d2o_OB,val_smooth);

    cl3_corr=f_ringCorr(cl3_corr,val_smooth);
    cl3_180_corr=f_ringCorr(cl3_180_corr,val_smooth);
    cl3_OB_corr=f_ringCorr(cl3_OB,val_smooth);

    empty_corr=f_ringCorr(empty_corr,val_smooth);
    empty_180_corr=f_ringCorr(empty_180_corr,val_smooth);
    empty_OB_corr=f_ringCorr(empty_OB,val_smooth);    
    
    
    
elseif ind_debug==0,


    d2o     =f_DC(d2o,DC);
    d2o_180 =f_DC(d2o_180,DC);
    d2o_OB  =f_DC(d2o_OB,DC);

    cl3      =f_DC(cl3,DC);
    cl3_180  =f_DC(cl3_180,DC);
    cl3_OB   =f_DC(cl3_OB,DC);

    empty      =f_DC(empty,DC);
    empty_180  =f_DC(empty_180,DC);
    empty_OB   =f_DC(empty_OB,DC);


    % flat field correction
    d2o =f_ffcorrect(d2o_OB ,d2o );
    d2o_180 =f_ffcorrect(d2o_OB ,d2o_180 );

    cl3 =f_ffcorrect(d2o_OB ,cl3 );
    cl3_180 =f_ffcorrect(cl3_OB ,cl3_180 );

    empty =f_ffcorrect(empty_OB ,empty );
    empty_180 =f_ffcorrect(empty_OB ,empty_180 );



    %dose correction
    d2o =f_doseCorrect(d2o ,xdose,ydose);
    d2o_180 =f_doseCorrect(d2o_180 ,xdose,ydose);
    d2o_OB =f_doseCorrect(d2o_OB,xdose,ydose);

    cl3 =f_doseCorrect(cl3 ,xdose,ydose);
    cl3_180 =f_doseCorrect(cl3_180 ,xdose,ydose);
    cl3_OB =f_doseCorrect(cl3_OB,xdose,ydose);

    empty =f_doseCorrect(empty ,xdose,ydose);
    empty_180 =f_doseCorrect(empty_180 ,xdose,ydose);
    empty_OB =f_doseCorrect(empty_OB,xdose,ydose);
    
        %% experimental: ring artefact correction
    d2o=f_ringCorr(d2o,val_smooth);
    d2o_180=f_ringCorr(d2o_180,val_smooth);
    d2o_OB=f_ringCorr(d2o_OB,val_smooth);

    cl3=f_ringCorr(cl3,val_smooth);
    cl3_180=f_ringCorr(cl3_180,val_smooth);
    cl3_OB=f_ringCorr(cl3_OB,val_smooth);

    empty=f_ringCorr(empty,val_smooth);
    empty_180=f_ringCorr(empty_180,val_smooth);
    empty_OB=f_ringCorr(empty_OB,val_smooth);  
    
else
    disp('bad ind_debug')
end

