    try
       a=rand(3)
       b=rand(5)
       c=[a,b]
    catch e %e is an MException struct
        fprintf(2,'%s\n',e.identifier);
        fprintf(2,'%s\n',e.message);
        % more error handling...
    end