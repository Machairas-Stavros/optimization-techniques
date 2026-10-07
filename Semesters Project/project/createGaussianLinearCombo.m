function gsns_res = createGaussianLinearCombo(u1_grid,u2_grid,m1,m2,s1,s2,c)

    % This function implements the calculation of the analytic expression
    % of the approximation surface.
    %
    % u1_grid and u2_grid is the result between the meshgrid implementation
    % between u1 and u2 vectors and represent the total number of
    % combinations between u1 and u2 points. Its size is N1xN2 where u1 has
    % N1 points and u2 N2 points.
    %
    % m1 is an NGx1xNSP matrix including the mean values across u1.
    % m2 is an NGx1xNSP matrix including the mean values across u2.
    %
    % s1 is an NGx1xNSP matrix including the sigma values across u1.
    % s2 is an NGx1xNSP matrix including the sigma values across u2.
    %
    % c is an NGx1xNSP matrix including the constant values for each
    % gaussian.
    %
    % ng is the total number of gaussians used.
    %
    % nsp is the number of samples that exist in the population.
    %
    % gsns_res is a 3D matrix of size N1xN2xNSP including the estimated
    % surface for each chromosome.

    ng = size(m1,1);
    
    gsns = {};
    
    for i=1:ng
        gsns{end+1} = c(i,:,:).*exp(-(u1_grid-m1(i,:,:)).^2./(2*s1(i,:,:).^2)-...
            (u2_grid-m2(i,:,:)).^2./(2*s2(i,:,:).^2));
    end
    
    gsns_res = zeros(size(u1_grid,1),size(u1_grid,2),size(m1,3));

    for i=1:ng
        gsns_res = gsns_res+gsns{i};
    end

    clear gsns
end
