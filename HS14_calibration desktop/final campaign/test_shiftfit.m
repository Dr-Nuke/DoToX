 figure(1);clf();
 plane=900;
 ref=recon(:,:,plane,1);
 im=-f.fraccircshift(recon2(:,:,plane,8),[-0.05,0.09]);
 im=recon2(:,:,plane,8);
 imshow((im)',[-0.005,0.005]);set(gca,'Ydir','normal')
 colorbar