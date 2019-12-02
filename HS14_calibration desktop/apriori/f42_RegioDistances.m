function d = f42_RegioDistances(t,regio,l)
% this function calculates the distance per region /material that a beam
% travelled through.
% xy = intersection coordinates
% regio = region swich array : 1:air (outside) 2:AL 3:D2 4: inside
% l = line object 

dist=([t;1]-[0;t])*l.l; %distances of segments

d=[0;0;0;0;0];
for i = 1:5
    d(i)=sum(dist(regio==i)); % cumulated distances
    
end

