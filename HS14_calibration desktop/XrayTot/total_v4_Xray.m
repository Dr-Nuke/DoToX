%% this is the 4th recon version, made for xrays
%% This script is the total X-ray script for the Zboray campaign 16.12.2016


%%
%init script
s4_init

%get name of thisscript
M.pwd=pwd; %get m file name

% manual settings
M.proofs=1;
M.d.rawfolder='C:/data/20161216 Campaign 5 Robert Z/';
M.d.nfiles=18;
M.d.imsize=[1024 640];
M.d.header_size=2048; 
M.d.filenames={ '11h28 Tmax 70kV 10mA.seq';...
            '11h31 Tmax 90kV 9mA.seq';...
            '11h34 Tmax 110kV 5mA.seq';...
            '12h04 B2 70kV 10mA.seq';...
            '12h06 B2 90kV 9mA.seq';...
            '12h08 B2 110kV 5mA.seq';...
            '12h33 B3 70kV 10mA.seq';...
            '12h34 B3 90kV 9mA.seq';...
            '12h35 B3 110kV 5mA.seq';...
            '12h59 B4 70kV 10mA.seq';...
            '13h01 B4 90kV 9mA.seq';...
            '13h03 B4 110kV 5mA.seq';...
            '13h37 leer 70kV 10mA.seq';...
            '13h39 leer 90kV 9mA.seq';...
            '13h41 leer 110kV 5mA.seq';...
            '14h14 voll 70kV 10mA.seq';...
            '14h16 voll 90kV 9mA.seq';...
            '14h18 voll 110kV 5mA.seq';...
            '20161220 DC 50 frames 30Hz.seq'};

cases=[1:18];
disp(sprintf('processing cases %s',num2str(cases)))
t(1)=toc;
for i=cases
disp(sprintf('case %02d\n',i))

M.i=i;
    
% load data       
M = f4_loadVideo(M,i);


% find all the relevant frames and dispose the rest
% manual guess
M.d.crop1=[  20,1354;    %1 % kepp only thoes frames (including)  
            25,1429;    %2
            24,1427;    
            29,1432;
            28,1431;    %5
            34,1438;
             4,1411;
            14,1418;
            16,1419;
            28,1431;    %10
            19,1424;
            16,1420;
            12,1416;
            14,1417;
            11,1414;    %15
            17,1419;
            14,1417;
            24,1427;
             1,50;];


M=f4_CropFrames(M,i); % CropFrames, returns M.im

        
       %% 
        
M.d.axlim=[ 70,130,400,440, %1 % required for the small axes
            70,130,400,440,
            70,130,400,440,
            70,130,400,440,
            70,130,400,440, %5
            70,130,400,440,
            70,130,400,440,
            70,130,400,440,
            70,130,400,440,
            70,130,400,440, %10
            70,130,400,440,
            70,130,400,440,
            70,130,400,440,
            70,130,400,440,
            70,130,400,440, %15
            70,130,400,440,
            70,130,400,440,
            70,130,400,440,
            70,130,400,440, ]; 
        

%% correct for the beam non-uniformity
% manual input of the dose areas
M.dose=[[5,427,429];...     % x start value
    [162,630,627];...       % x end value
    [100,100,700];...       % y start value
    [1010,229,1014]];       % y end value

% make the dose areas and check them
M=f4_DoseArea(M);
 
%% calculate & apply the beam non uniformity
M=f4_FitBeam(M);

% check upon it
f4_FitBeamCheck(M);
M.sinoslice{i}=squeeze(imc(:,45,:));
%%
M.imcname{i}=sprintf('imc%02d',i);
eval([M.imcname{i} '= M.imc;']);
fprintf('\n saving %s...',M.imcname{i})
savefast(strcat(M.pwd,'\',M.imcname{i}),M.imcname{i});
eval(['clear ' M.imcname{i}]); % delete from workspace
fprintf('done. \n')
%%

fprintf('freeing up memory...')
M=rmfield(M,'unibeam');
M=rmfield(M,'raw');
M=rmfield(M,'im');
M=rmfield(M,'imc');
 M=rmfield(M,'f'); % more for iteration reasons
fprintf('done. \n')

t(i+1)=toc;
dt=t(i+1)-t(i);
disp(sprintf('case %02d elapsed time: %5.f sec/ %5.2f min /%4.2f h',i,dt,dt/60,dt/3600))

end
save('M','M')
disp(fprintf('total time %5. min',toc/60))
%%
% reference pairs
M.refpairs=[2,14;
            3,15;
            4,13;
            5,14;
            6,15;
            7,13;
            8,14;
            9,15;
            10,13;
            11,14;
            12,15;
            16,13;
            17,14;
            18,15;];
          
        
% matching 
for i =1:size(M.refpairs,1)
    M.i=i;
    M=f4_SinoMatch(M);
end
% check matches
for i =1:size(M.refpairs,1)
    M.i=i;
    f4_SinoMatchCheck(M);
end



%% go towards reconstruction
% cut away undesired stuff
M.rpx=[145,435,50,1015]; % reconpix x1 x2 y1 y2
M.recsize=291;
for i =M.refpairs(:,1)
    M.i=i
    M = f4_centering(M);


end
%coupling to old LFT script
reconblock=M.recon;
save('reconblock','reconblock')
%%
s3_LFT2






