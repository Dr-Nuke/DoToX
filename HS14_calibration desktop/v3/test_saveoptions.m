% tic
% load('/media/data/cbolesch/Tomo_HS14/processed/block_2blc_d2o_.mat')
% toc

tic
a=rand(2450,1200,400,'single');
t2=toc;
disp(strcat('creating random:_ ',num2str(t2)));

imin=[1,101,201,301];
imax=[100,200,300,400];
subblock=zeros(2450,1200,100,'single');


for i=1:4
    t1=toc
    filename=strcat('testblock.mat',num2str(i));
    filepath=strcat('/media/data/cbolesch/Tomo_HS14/processed/',filename);
    subblock=a(:,:,imin(i):imax(i));
    save(filepath,'subblock','-v7.3');
    t2=toc;
    disp(strcat('saving_',filename,'_',num2str(t2-t1)));
end

tic
save('/media/data/cbolesch/Tomo_HS14/processed/testblock_full.mat','a','-v7.3');
toc