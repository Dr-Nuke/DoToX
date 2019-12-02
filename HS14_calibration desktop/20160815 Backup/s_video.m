% % create a video vom image files
% for i=1:size(d2o,3)
%     f(i)=im2frame(d2o_180/max(d2o_180(:)),gray);%d2o/(max(d2o(:)-min(d2o(:)))-min(d2o(:))));
% end
% movie(f)

% matrix2avi(d2o*5,'file','C:\data\tomo_HS14\processed\2.avi','map',gray)

mov(1:375) = struct('cdata', [],'colormap', []);


for i=1:375
   mov(1,i).cdata=d2o_dose1(:,:,i)*100;
    
end
movie2avi(mov, 'C:\data\tomo_HS14\processed\RGB3.avi', 'compression', 'None');