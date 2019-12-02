function [M] = f4_DMT(M)
%% this function sets up the meta data table for the zboray xray campaign
xlspath='G:/cbolesch/Data/';
xlsfname='Data Overview.xlsx';
xlsfile=strcat(xlspath,xlsfname);

formatSpec = '%{yyyy-MM-dd}D%{HH:mm}D%s%s%s%s%s%s%f%f%f%f%s';
DMT=readtable(xlsfile,'Range','A3:N75');
DMT(1:53,:)=[];     %discard non interesting stuff
DMT(:,{'facility','Var4','tomo_radio','rawDataPath','videoPath'})=[];


%fileNumbers=
fileNames={ '11h28 Tmax 70kV 10mA.seq';...
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
FileNumber=[1:19]';        

fnames=table(fileNames,'VariableNames',{'File_Names'});

DMT=[fnames,DMT];

M.DMT=DMT;
end

