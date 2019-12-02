%s3_filter tester



addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop')
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\v2_1')
load('C:\data\Tomo_HS14\processed\block_2blc_dc__1.mat')
load('C:\data\Tomo_HS14\processed\block_2blc_cl3_4.mat')

a=block_2blc_dc__1(:,:,1);
%imbo3(a,1);

%% line filter
for i=1
%a=block_2blc_dc__1(:,:,3);
%a=block_2blc_cl3_4(:,:,1);
a=squeeze(block_2blc_cl3_4(:,:,1));
thresh_medi1=0.05; % for gamma spots (single pixel)
thresh_medi2=0.15; % for large spots DC 0.2 good value
kernel_medi2=f_KernelGen(9,9,7.5);
%b=f3_MediFilter2(a,kernel_medi2,thresh_medi2,'median' );


% filter gammaspots(median)
% after-filtering to remove circles left overy by spotfilter
% b=f_MediFilter1(a,thresh_medi1);    %vertical
%b=f_MediFilter1_2(b,thresh_medi1);  % horizontal (or vice versa)
b=a;%f3_CamLineFiltDC(a);


start=1;
stop=2^16;

bins=linspace(start,stop,1000);


figure(1)
cla
ha=histogram(a(:),bins,'EdgeColor','none');
ha.FaceColor = [1 0 0];
hold on
hb=histogram(b(:),bins,'EdgeColor','none','facealpha',0.5);
legend('original','filtered')
hb.FaceColor = [0 1 0];
set(gca, 'YScale', 'log')
diff=ha.Values-hb.Values;
f=sum(abs(diff))/prod(size(a));
text((stop-start)/2,max(ha.Values)/2,strcat(num2str(f*10),' percent of pixels modified'))
grid on
%plot(bins(1:end-1),diff,'r')

aa=log(a);
bb=log(b);
aamin=9.5;%min(bb(:));
aamax=max(bb(:));




c=a==b;
figure(2)
ax1=subplot(1,3,1);
cla
imagesc(squeeze(aa)');colormap(jet);axis equal; axis tight;set(gca,'YDir','normal');

caxis([aamin aamax])
%colorbar
ax2=subplot(1,3,2);
cla
imagesc(squeeze(bb)');colormap(jet);axis equal; axis tight;set(gca,'YDir','normal');
caxis([aamin aamax])
%colorbar
ax3=subplot(1,3,3);
imagesc(squeeze(c)');colormap(jet);axis equal; axis tight;set(gca,'YDir','normal');

linkaxes([ax1,ax2,ax3],'xy')
% truesize

end



%% spot filter
for i =1 
%a=block_2blc_dc__1(100:end,:,3);
%a=block_2blc_cl3_4(:,:,1);
a=squeeze(cl3(:,10,:));
thresh_medi1=0.05; % for gamma spots (single pixel)
thresh_medi2=0.15; % for large spots DC 0.2 good value
kernel_medi2=f_KernelGen(9,9,7.5);
%b=f3_MediFilter2(a,kernel_medi2,thresh_medi2,'median' );


% filter gammaspots(median)
% after-filtering to remove circles left overy by spotfilter
b=f_MediFilter1(a,thresh_medi1);    %vertical
%b=f_MediFilter1_2(b,thresh_medi1);  % horizontal (or vice versa)



start=1;
stop=2^16;

bins=linspace(start,stop,1000);


figure(1)
cla
ha=histogram(a(:),bins,'EdgeColor','none');
ha.FaceColor = [1 0 0];
hold on
hb=histogram(b(:),bins,'EdgeColor','none','facealpha',0.5);
legend('original','filtered')
hb.FaceColor = [0 1 0];
set(gca, 'YScale', 'log')
diff=ha.Values-hb.Values;
f=sum(abs(diff))/prod(size(a));
text((stop-start)/2,max(ha.Values)/2,strcat(num2str(f*10),' percent of pixels modified'))
grid on
%plot(bins(1:end-1),diff,'r')

aa=log(a);
bb=log(b);
aamin=9.5;%min(bb(:));
aamax=max(bb(:));
c=a==b;
figure(2)
ax1=subplot(3,1,1);
cla
imagesc(squeeze(aa)');colormap(jet);axis equal; axis tight;set(gca,'YDir','normal');

caxis([aamin aamax])
%colorbar
ax2=subplot(3,1,2);
cla
imagesc(squeeze(bb)');colormap(jet);axis equal; axis tight;set(gca,'YDir','normal');
caxis([aamin aamax])
%colorbar
ax3=subplot(3,1,3);
imagesc(squeeze(c)');colormap(jet);axis equal; axis tight;set(gca,'YDir','normal');

linkaxes([ax1,ax2,ax3],'xy')
% truesize


end

%% DC spot filter
for i =1 
a=block_2blc_dc__1;
b=a;
thresh_medi3=0.05; % for gamma spots (single pixel)
kernel=[1,5]
[x,y,z]=size(a);



    for j=1:y % iterate y-planes
        
        b(:,j,:)=f3_MediFilter3(squeeze(a(:,j,:),thresh_medi3,)
        
        
    
    end


    

b=f_MediFilter1(a,thresh_medi1);    %vertical




start=1;
stop=2^16;

bins=linspace(start,stop,1000);


figure(1)
cla
ha=histogram(a(:),bins,'EdgeColor','none');
ha.FaceColor = [1 0 0];
hold on
hb=histogram(b(:),bins,'EdgeColor','none','facealpha',0.5);
legend('original','filtered')
hb.FaceColor = [0 1 0];
set(gca, 'YScale', 'log')
diff=ha.Values-hb.Values;
f=sum(abs(diff))/prod(size(a));
text((stop-start)/2,max(ha.Values)/2,strcat(num2str(f*10),' percent of pixels modified'))
grid on
%plot(bins(1:end-1),diff,'r')

aa=log(a);
bb=log(b);
aamin=9.5;%min(bb(:));
aamax=max(bb(:));
c=a==b;
figure(2)
ax1=subplot(3,1,1);
cla
imagesc(squeeze(aa)');colormap(jet);axis equal; axis tight;set(gca,'YDir','normal');

caxis([aamin aamax])
%colorbar
ax2=subplot(3,1,2);
cla
imagesc(squeeze(bb)');colormap(jet);axis equal; axis tight;set(gca,'YDir','normal');
caxis([aamin aamax])
%colorbar
ax3=subplot(3,1,3);
imagesc(squeeze(c)');colormap(jet);axis equal; axis tight;set(gca,'YDir','normal');

linkaxes([ax1,ax2,ax3],'xy')
% truesize


end


%