figure(76)
clf
hold on
scrange=380:2200
for i=1:4
    scatter3(c(scrange,i,1),c(scrange,i,2),scrange,'.');
end
set(gca, 'CameraViewAngle', 8.4)