%% calculate the reconstruction data
if ind_debug==1,

    raw_cl3_d2o=-log(cl3_corr./d2o_corr);
    raw_cl3_d2o_180=-log(cl3_180_corr./d2o_180_corr);

%     raw_cl3_emp=-log(cl3_corr./empty_corr);
%     raw_cl3_emp_180=-log(cl3_180_corr./empty_180_corr);

    raw_d2o_emp=-log(d2o_corr./empty_corr);
    raw_d2o_emp_180=-log(d2o_180_corr./empty_180_corr);

    % fix zeros
elseif ind_debug==0,
    
    cl3=-log(cl3 ./d2o );
    cl3_180=-log(cl3_180 ./d2o_180 );

    d2o=-log(d2o ./empty );
    d2o_180=-log(d2o_180 ./empty_180 );
else
    disp('bad ind_debug')
end