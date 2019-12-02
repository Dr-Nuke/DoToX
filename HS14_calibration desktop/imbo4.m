function imbo4( im)

imshow(squeeze(im)',[]);
colormap(gray); 
axis equal; 
axis tight;
set(gca,'YDir','normal'); 
colorbar
end

