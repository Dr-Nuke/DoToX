% test the astra toolbox.

% first load some data
 load('T')
 d=f.BoLoad(T.fnames.add{6},T);
% 
 dat=d(173:491,:,:);
%%
plane=640;
sinogram=-log(squeeze(dat(:,plane,:)));
sinogram=f.fraccircshift(sinogram,-T.Cen.fitshift(2,plane))';


tic;
for i=1:10
rec2=a.FBP(sinogram,A);
end
toc
