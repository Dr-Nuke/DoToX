figure(1);clf; hold on;

a=0;                            % start of interval
b=1;                            % end of interval
epsilon=0.001;               % accuracy value
iter= 10;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k=0;                            % number of iterations

g=@(x)-(1-(x-0.78).^2);             % quest function
x1n=@(a,b,tau)a+(1-tau)*(b-a)   % lefz mid point
x2n=@(a,b,tau)a+tau*(b-a);      % riht mid point

x1=x1n(a,b,tau);             % computing x values
x2=x2n(a,b,tau);
xlist=[x1,x2];

f_x1=g(x1);                     % computing values in x points
f_x2=g(x2);
flist=[f_x1,f_x2];

plot(x1,f_x1,'kx')              % plotting x
plot(x2,f_x2,'kx')
xlim([a,b])
grid on
col=jet(iter)
while ((abs(b-a)>epsilon) && (k<iter))
    k=k+1
    if(f_x1>f_x2)
        b=x2;
        x2=x1;
        x1=a+(1-tau)*(b-a);
        xlist=[xlist,x1];
        f_x2=f_x1;
        
        
        f_x1=g(x1);
        flist=[flist,f_x1];
        plot(x1,f_x1,'x','color',col(k,:));

    else
        a=x1;
        x1=x2;
        x2=a+tau*(b-a);
        xlist=[xlist,x2];
        
        f_x1=f_x2;
        f_x2=g(x2);
        flist=[flist,f_x2];
        
        plot(x2,f_x2,'+','color',col(k,:))
    end
    
    
end


% chooses minimum point
if(f_x1<f_x2)
    sprintf('x_min=%f', x1)
    sprintf('f(x_min)=%f ', f_x1)
    plot(x1,f_x1,'ro')
else
    sprintf('x_min=%f', x2)
    sprintf('f(x_min)=%f ', f_x2)
    plot(x2,f_x2,'ro')
end
title('')