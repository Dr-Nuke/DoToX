function M = f4_SinoMatch(M)
i=M.i;
ii=M.refpairs(i,1);


sino=M.sinoslice{ii};
ref=M.sinoslice{M.refpairs(i,2)};


[sino2,ref2]=f4_SinoPad(sino,ref);


%% iterate 2d shift grid

sxb=linspace(-5,5,21);  % basis iteration grd
sfb=linspace(-5,5,21);  % i is i-dimenion, f is frames-dimension
sx=sxb;                 % actual iteration grid
sf=sfb; 
h=figure(10);
iterations=3;
for j=1:iterations % number of iterations

    diff2d=zeros(length(sx),length(sf)); % preallocate map
    for k=1:length(sx) % x dimension
        for l=1:length(sf) % frame dimension
            m=l+(k-1)*length(sx);
            f_BoCount(m,20,10,5)
            im=f4_ShiftXY(sino2,[sx(k),sf(l)]);
            diff2d(k,l)=sum(abs(ref2(:)-im(:)));
        end
    end
    fprintf('\n')
    % optimize for region
    [Min,I] = min(diff2d(:)); % find shift & index of minimum difference 
    [I_row, I_col] = ind2sub(size(diff2d),I); % convert index to 2d index
    
    % plot
    ax(j)=subplot(1,iterations,j);
    imagesc(diff2d','XData', sx, 'YData', sf);set(gca,'YDir','normal');
    colorbar; colormap gray
    xlabel('x')
    if j==1
        ylabel('y')
    end
    
    % make title for plot
    titlestring=sprintf('iter. %d, min dev. %.1f \n shift: %4.2f %4.2f',...
        j,diff2d(I_row,I_col),sx(I_row),sf(I_col));
    if j==2
        titlestring=sprintf('case %d \n %s',ii,titlestring);
    end
    title(titlestring)
    axis equal
    axis tight
    
    if j==3 % assign the shift to M
        M.shift(ii,:)=[sx(I_row),sf(I_col)];
    else % zoom in the shifts
        sxb=sxb/10;
        sfb=sfb/10;
        sx=sxb+sx(I_row);
        sf=sfb+sf(I_col);
    end
end
    fname=sprintf('Fig10 %02d sino shift match',ii);
    f_BoFig2PDF(h,fname) 
%M.shift(i,:)=[]
end

