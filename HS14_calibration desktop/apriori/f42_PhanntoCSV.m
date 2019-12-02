%% Thsi script writes the phantom to a csv file

P=f42_ChannelPhanTransform;

%make a table out of the struct


Varnames={'Nr','C_L_swich','center_XY','Radius','Start_Angle',...
    'Angular_Distance','Region_swich','Line_start_XY','Line_end_XY'};

T=table(P.n',P.CL',P.c,P.r',P.p0',P.dp',P.reg,P.Lx1,P.Lx2,...
    'VariableNames',Varnames)

C=T(T{:,2}==1,:);% only cirles
C.Nr=(1:length(C.Nr))';
C=C(:,[1,3,4,5,6,7]);

L=T(T{:,2}==0,:);%only lines
L.Nr=(1:length(L.Nr))';
L=L(:,[1,7,8,9]);
writetable(C,'Phantom_Circles.csv')
writetable(L,'Phantom_Lines.csv')