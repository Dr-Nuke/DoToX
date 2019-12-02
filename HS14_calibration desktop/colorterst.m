figure()
hold on
a=hsv;
for i=1:100
    ind_col=ceil(size(a,1)*i/100);
    plot([1,2],[3+i/10,4],'color',a(ind_col,:))
    for j=1:5
       pause(0)
    end
     fprintf('%d ',j)
end
 
;