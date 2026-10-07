function func = evaluationFunct(u1_grid,u2_grid)

    % This function calculates the surface of the evaluation function.
    %
    % u1_grid and u2_grid is the grid of values containing all the possible
    % points combinations between the vectors u1 and u2. It derives from
    % the meshgrid function.
    %
    % func is the surface of the evaluation function.
    
    func = sin(u1_grid+u2_grid).*sin(u2_grid.^2);

end