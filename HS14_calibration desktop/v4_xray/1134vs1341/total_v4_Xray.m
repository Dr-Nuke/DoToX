%% this is the 4th recon version, made for xrays

%get name of thisscript
[~,name,~]=fileparts(mfilename('fullpath')); %get m file name

%init script
s4_init

% flag for displaying the proofs
M.proofs=1;

% load film movie
load('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v4_xray/1134h/imc2.mat')

% load empty movie
load('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v4_xray/1341h/imc1.mat')
load('match')
%imc1=double(imc1);

imc2=imc2(:,:,1:1455); %crop
%imc2=double(imc2);


%% shift the stack
imc2c=zeros(size(imc2));
for i=1:size(imc2,2)
    f_BoCount(i,20,10,5)
    imc2c(:,i,:)=reshape(f4_ShiftXY(squeeze(imc2(:,i,:)),match),[640,1,1455]);
end
%%
imc2c=single(imc2c);
M.raw=-log(imc2c./imc1);

M.proofs=1;


%% find all the relevant frames and dispose the rest

% manual guess
M.d.crop=[30,1386]; % kepp only thoes frames (including)
M.d.axlim=[70,130,400,440]; % required for the small axes

M=f4_CropFrames(M); % CropFrames, returns M.im

%% correct for the beam non-uniformity
% manual input of the dose areas
M.dose=[[5,427,429];...     % x start value
    [162,630,627];...       % x end value
    [100,100,700];...       % y start value
    [1010,229,1014]];       % y end value

% % make the dose areas and check them
% M=f4_DoseArea(M);
% 
% %% calculate & apply the beam non uniformity
% M=f4_FitBeam(M);
% 
% % check upon it
% f4_FitBeamCheck(M);


%% go towards reconstruction
% cut away undesired stuff

% make circshift ~60deg
M.rpx=[145,435,50,1015]; % reconpix x1 x2 y1 y2
M.recsize=291;

M = f4_centering(M);

%coupling to old LFT script
reconblock=M.recon;
save('reconblock','reconblock')
%%
s3_LFT2






