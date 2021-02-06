t = 0:.1:1000;
y = -1000:.1:1000;
Z = 0:.1:1000; 

ismon = true; 
for t = .1:.1:1000
    for y = 1:1:1000
        fun = @(Z) exp(-Z.^4).*sqrt(Z*y).*besselj(nu, Z*y/t);
        if exist('curr', 'var') == false
            curr = integral(fun, 0, 1000); 
            continue;
        end
        if exist('nextcurr', 'var') == false
            nextint = integral(fun, 0, 1000); 
            if nextint - curr >=0
                mon = 1;
            else
                mon = 2; 
            end
            continue; 
        end
        nextint = integral(fun, 0, 1000); 
        if nextint - curr < 0 && mon == 1
            keyboard;
        end
        if nextint - curr >= 0 && mon == 2
            keyboard;
        end

        
    end
end
