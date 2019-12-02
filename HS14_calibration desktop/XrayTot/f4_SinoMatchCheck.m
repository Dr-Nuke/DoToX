function []= f4_SinoMatchCheck(M)
%checks upon the Sino match
i=M.i;
ii=M.refpairs(i,1);

xwin=[1,20;
    240,280;
    520,560];

ywin=[435,454;
        150,190;
        430,470];



sino=M.sinoslice{M.refpairs(i,1)};
ref=M.sinoslice{M.refpairs(i,2)};


    h=figure(11);clf;
for j=1:3
    
    l=(j-1)*3+1;

    ax(l)=subplot(3,3,l);
    imshowpair(sino,ref)
    xlim(xwin(j,:))
    ylim(ywin(j,:))
    if j==1
        title(sprintf('unshifted'))
    end

    ax(l+1)=subplot(3,3,l+1);
    imshowpair(f4_ShiftXY(sino,M.shift(ii,:)),ref)
    xlim(xwin(j,:))
    ylim(ywin(j,:))
    if j==1
        title(sprintf('case %d %4.2f %4.2f \n XY shift',...
            ii,M.shift(ii,1),M.shift(ii,2)))
    end

    ax(l+2)=subplot(3,3,l+2);
    imshowpair(f4_ShiftXY(sino,[0,M.shift(ii,2)]),ref)
    xlim(xwin(j,:))
    ylim(ywin(j,:))
    if j==1
        title(sprintf('only Y shift'))
    end
end
    fname=sprintf('Fig11 %02d sino smatch check',ii);
    f_BoFig2PDF(h,fname) 

end

