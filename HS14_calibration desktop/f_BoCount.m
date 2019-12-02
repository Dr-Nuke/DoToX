function [ output_args ] = f_BoCount(c,f,r,d)
% puts counter display values in order
% f_BoCount(c,f,r,d)
% c = counter
% f = only every f-th number shown 
% r = output line width in numbers
% d = digits per number

if mod(c,f)==0;
    string=strcat('%',num2str(d),'.d');
    if mod(c,r*f)==0
       fprintf(strcat(string,'\n'),c)
        
    else
       fprintf(string,c)
    end
end

end

