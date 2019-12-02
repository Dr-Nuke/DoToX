p00 = zeros(8,9,1000);
p10 = zeros(8,9,1000);
p01 = zeros(8,9,1000);
p20 = zeros(8,9,1000);
p11 = zeros(8,9,1000);
p02 = zeros(8,9,1000);

%%
figure(1);clf
ax(1)=gca;
hold on
grid on
title('p00')

figure(2);clf
ax(2)=gca;
hold on
grid on
title('p01')

figure(3);clf
ax(3)=gca;
hold on
grid on
title('p10')

figure(4);clf
ax(4)=gca;
hold on
grid on
title('p20')

figure(5);clf
ax(5)=gca;
hold on
grid on
title('p11')

figure(6);clf
ax(6)=gca;
hold on
grid on
title('p02')

for cas=[1,6,7,8]
    disp(cas)
    for rep=1:4
        for frame=80:(80+1357);
            p00(cas,rep,frame) = T.fit.fit{cas,rep}{frame}.p00;
            p10(cas,rep,frame) = T.fit.fit{cas,rep}{frame}.p10;
            p01(cas,rep,frame) = T.fit.fit{cas,rep}{frame}.p01;
            p20(cas,rep,frame) = T.fit.fit{cas,rep}{frame}.p20;
            p11(cas,rep,frame) = T.fit.fit{cas,rep}{frame}.p11;
            p02(cas,rep,frame) = T.fit.fit{cas,rep}{frame}.p02;
        end
    end
    plot(ax(1),squeeze(p00(cas,rep,:)),'Displayname',num2str(cas));
    xlim([80,(80+1357)])
    plot(ax(2),squeeze(p10(cas,rep,:)),'Displayname',num2str(cas));
    xlim([80,(80+1357)])
    plot(ax(3),squeeze(p01(cas,rep,:)),'Displayname',num2str(cas));
    xlim([80,(80+1357)])
    plot(ax(4),squeeze(p20(cas,rep,:)),'Displayname',num2str(cas));
    xlim([80,(80+1357)])
    plot(ax(5),squeeze(p11(cas,rep,:)),'Displayname',num2str(cas));
    xlim([80,(80+1357)])
    plot(ax(6),squeeze(p02(cas,rep,:)),'Displayname',num2str(cas));
    xlim([80,(80+1357)])
end

for i=1:6
    ax(i);
    legend();
end
%%
masterfit=T.fit.fit{1,1}{1}
masterfit.p00=mean(p00(p00~=0));
masterfit.p01=mean(p01(p01~=0));
masterfit.p10=mean(p10(p10~=0));
masterfit.p20=mean(p20(p20~=0));
masterfit.p02=mean(p02(p02~=0));
masterfit.p11=mean(p11(p11~=0))

[x,y]=ndgrid(1:640,1:1024);
masterbeam=feval(masterfit,x,y);
testbeam=feval(T.fit.fit{1,1}{1},x,y);
masterbeam=masterbeam/mean(masterbeam(:));
save('masterbeam','masterbeam')








