function [sino,ref] = f4_SinoPad(sino,ref)
% returns the smaller sinogram padded wit hthe last fram 
% such that both have same size afterwards

s=size(sino);
r=size(ref);
if s(1)~=r(1)
    error('two sinograms are not equally sized in x direction')
end

if s(2)<r(2) 
    sino=[sino,repmat(sino(:,end),[1,r(2)-s(2)])];
elseif r(2)<s(2) 
    ref=[ref,repmat(ref(:,end),[1,s(2)-r(2)])];
end




end

