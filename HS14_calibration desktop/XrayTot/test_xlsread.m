xlspath='G:/cbolesch/Data/';
xlsfname='Data Overview.xlsx';
xlsfile=strcat(xlspath,xlsfname);

formatSpec = '%{yyyy-MM-dd}D%{HH:mm}D%s%s%s%s%s%s%f%f%f%f%s';
DMT=readtable(xlsfile,'Range','A3:N75');
DMT(1:53,:)=[];     %discard non interesting stuff
DMT(:,{'facility','Var4','tomo_radio','rawDataPath','videoPath'})=[];

% for i =4:8
%     a=zeros(1,19);
%     for j=1:19
%         a(j)=str2num(DMT{j,i}{1});
%     end
%     DMT(:,10)=table(a');
%     DMT(:,4)=[]
% end
    
DMT