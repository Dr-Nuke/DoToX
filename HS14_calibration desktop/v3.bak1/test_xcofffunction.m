% test for the new function to encapsulate the xcov subgrid  shift
clear a b 

testStpIdx=[2,4,4] % need the block2 of emp and blocks 4 for cl3 and d2o

for j= 1


%     block_old=strcat(nam_block,nam_stp{testStpIdx(j)},'_',pfx_cas{j});
%     path_block=strcat(dir_work,block_old);
%     tiltblock=f3_loadSingleVariableMATFile(path_block);



    for i=1:374
        a(j,i)=f3_subpixcorr(tiltblock(:,600,i),tiltblock(:,2000,i));
    end

end

figure(6);clf;
for j=1
    plot(a(j,:)); hold on
end

%hold on;
%plot(b);

xlabel('projecion number')
ylabel('tilt in pixel')
%legend({'subpixel shift','integer pixel tilt'})
legend({'empty','d20','cl3'})
meantilt2=mean(a(:))
title(strcat('tilt in pixel. mean tilt: ',' ', num2str(meantilt2)))
xlim([1,374])

