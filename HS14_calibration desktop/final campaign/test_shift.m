% we have d, d_add1
cas=8;rep=1;T.cas=cas;T.rep=rep;
dn=f.BoLoad(T.fnames.addNewshift{8},T); %new %w x-shift
do=f.BoLoad(T.fnames.add{8},T);     % old w/o x-shift
dr=f.BoLoad(T.fnames.addNewshift{1},T); %ref
%%


yplane=500;
sn=squeeze(dn(:,yplane,:));
so=squeeze(do(:,yplane,:));
sr=squeeze(dr(:,yplane,:));

crange=([-0.02,0.05])

h(1)=figure(1),clf;
set(h(1),'Position',[0,39,600,958]);
ax(1)=imshow(-log(sn./sr)',[crange]);set(gca,'YDir','normal');
title('new')

h(2)=figure(2),clf;
set(h(2),'Position',[518,40,600,958]);
ax(2)=imshow(-log(f.fraccircshift(sn,[-0.2,0])./sr)',[crange]);set(gca,'YDir','normal')
title('old')

h(3)=figure(3),clf;
set(h(3),'Position',[1032,40,600,958]);
ax(3)=imshow((sn-so)',[crange]);set(gca,'YDir','normal')
title('difference')

linkaxes(ax)

%% inside FindCentering

% lateral shift

for frame=1:(k180-5)
    [c,lags]=xcov(d(:,plane,1+frame),flipud(d(:,plane,k180+frame)),50,'coeff');
        if any(isnan(c)) %debug
            fprintf('%d, ',plane)
            continue
        else
            %disp(sprintf('%4d does not contain NaN',j))
            %ind_nan(k)=1;
        end
        [~,ind]=max(c);
        try
        ls(frame)=(lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/...
            (log(c(ind-1))-2*log(c(ind))+log(c(ind+1))))/2;
        catch
            ls(plane)=nan;
        end
end



%% 




