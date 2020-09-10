function pulse = get_pulse(N, per_nb)
    
    width = floor(N/per_nb/2);
    
    iteration = [zeros(width,1); ones(width, 1)];
    
    if width == N/per_nb
        pulse = kron(ones(per_nb, 1), iteration);
        return; 
    else
        pulse = kron(ones(per_nb, 1), iteration);
    end
    
    if floor(per_nb / 2) == 0
        pulse = [pulse; zeros(N - per_nb*width*2,1)];
    else
        pulse = [pulse; ones(N - per_nb*width*2,1)];
    end
        
end