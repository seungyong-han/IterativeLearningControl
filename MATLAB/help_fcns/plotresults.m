function p = plotresults(error_history, iteration_number, Algorithm, logscale, xshift, yshift, color)
    
    switch Algorithm
        case('RIA')
            name = 'Right Inverse Model Algorithm';
        case('LIA')
            name = 'Left Inverse Model Algorithm';
        case('IA')
            name = 'Inverse Model Algorithm';
        case('SDA')
            name = 'Steepest Descent Algorithm';
        case('SDA_ES')
            name = 'Suppression of the Eigenvalues';
        case('FD')
            name = 'Feedback Design Algorithm';
        otherwise
            name = 'none'
            fprintf(2, ['There is no algorithm ', Algorithm]);
    end

    
hold on
    if logscale
        set(gca, 'YScale', 'log')
    end
    
    nb = length(error_history);
    iteration_number = nb; 
    p = plot(1:nb, error_history, 'Color', color, 'DisplayName', name);
    plot([iteration_number, iteration_number], [0, error_history(nb)], '--', 'color', 'black');
    plot([0, iteration_number],[error_history(nb), error_history(nb)], '--', 'color', 'black');
    
    %text(iteration_number - 1000,error_history(nb)+1,['\textbf{', name, '}', char(10), 'Iteration number: ',num2str(iteration_number), char(10), 'error: ', num2str(error_history(nb))], 'interpreter', 'latex');
    text(iteration_number + xshift,error_history(nb)+yshift,['\textbf{', name, '}', char(10), 'Iteration number: ',num2str(iteration_number), char(10), 'error: ', num2str(error_history(nb))], 'interpreter', 'latex');
    ylabel('abs error');
    xlabel('iteration');
    
    plot(iteration_number, error_history(end), 'o', 'LineWidth', 8, 'Color', color);
    
hold off
end