function M = f4_FitBeam(M)
% fits the beam nonuniformity

mask=M.mask;
maskmesh=size(mask);

[x,y]=ndgrid(1:maskmesh(1),1:maskmesh(2));
xx=x;
yy=y;
xx(mask==0)=[];
yy(mask==0)=[];

if isfield(M,'imc')
    M=rmfield(M,'imc');
end

if isfield(M,'unibeam')
    M=rmfield(M,'unibeam');
end

M.unibeam=ones(size(M.im),'single');


% find the fits
disp('finding beam non-uniformity fits...')
for i=1:M.nFrame
    
    f_BoCount(i,50,10,5)
    im=squeeze(M.im(:,:,i));    % check out a slice from stack
    im(mask==0)=[];             % remove non-flatfield-areas
    M.f{i}= fit( [xx',yy'],double(im'), 'poly22' );   % make a fit
    M.unibeam(:,:,i)=feval(M.f{i},x,y);         % make the beam
    
end

% correct the non uniformity
M.imc=single(M.im)./M.unibeam;


end