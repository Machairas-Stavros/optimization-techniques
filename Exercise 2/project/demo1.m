clf
close all
clear
clc

rng(123)

syms x1 x2 gk  
varsvect = [x1;x2;gk];

f = x1^3*exp(-x1^2-x2^4);
f_handle = matlabFunction(f,'Vars',{'x1','x2'});
f_label = "f=x_1^3*e^{-x_1^2-x_2^4}";
x_values = -2.5:0.01:2.5;


[x1_grid_values,x2_grid_values] = meshgrid(x_values,x_values);
y_values_f = f_handle(x1_grid_values,x2_grid_values);

[minf,minfind] = min(y_values_f(:));
[x1_min_ind,x2_min_ind] = ind2sub(size(y_values_f),minfind);
x1_min = x1_grid_values(x1_min_ind,x2_min_ind);
x2_min = x2_grid_values(x1_min_ind,x2_min_ind);

figure()
surf(x1_grid_values,x2_grid_values,y_values_f)
title('Surface - '+f_label)
xlabel('x_1 Values')
ylabel('x_2 Values')
zlabel('f Values')
colorbar

figure()
plot3(x1_grid_values,x2_grid_values,y_values_f)
title('Plot3 - '+f_label)
xlabel('x_1 Values')
ylabel('x_2 Values')
zlabel('f Values')

figure()
contour(y_values_f,100)
title('Contour - '+f_label)
subtitle('100 Isolines')
xlabel('x_1 Values')
ylabel('x_2 Values')
colorbar

starting_points = [0,-1,1;0,-1,1];
gama_test_values = 0.05:0.05:1;

epsilon_sda = 10^-5;

fprintf('Execution of Steepest Descent Algorithm - Default g_k\n')
dg_sda = {};

for i=1:size(starting_points,2)
    dg_sda{end+1} = struct('StartingPoint',starting_points(:,i));
    mincalc = [];
    noreps = [];
    for j=1:length(gama_test_values)
        fprintf('Calculation of Minimum - StartingPoint=[%s] - g_k=%.2f\n',strjoin(string(starting_points(:,i)),','),gama_test_values(j))
        [mincalc(end+1,:),noreps(end+1,:)] = SteepestDescentAlgorithm(dg_sda{i}.StartingPoint,epsilon_sda,f,varsvect,gama_test_values(j));
    end
    dg_sda{i}.xValuesMin = mincalc;
    dg_sda{i}.NoReps = noreps;
end

for i=2:length(dg_sda)
    figure()
    subplot(2,2,[1 3])
    hold on
    for j=1:size(dg_sda{i}.xValuesMin,1)
        plot(dg_sda{i}.xValuesMin(j,1),dg_sda{i}.xValuesMin(j,2),'o')
    end
    hold off
    title('Variance of (x_1,x_2) Values in Function of g_k')
    xlabel('x_1 Values')
    ylabel('x_2 Values')
    legend("gk = "+num2cell(string(gama_test_values)))

    f_temp_calc = [];
    for j=1:size(dg_sda{i}.xValuesMin,1)
        funmininparg = num2cell(transpose(dg_sda{i}.xValuesMin(j,:)));
        f_temp_calc(end+1) = feval(f_handle,funmininparg{:});
    end

    subplot(2,2,2)
    plot(gama_test_values,f_temp_calc,'Marker','o','Color','Red')
    title('f(x_1,x_2) in Function of g_k')
    xlabel('g_k')
    ylabel('f(x_1,x_2)')

    subplot(2,2,4)
    plot(gama_test_values,dg_sda{i}.NoReps,'Marker','o','Color','Blue')
    title('Algorithm''s Repetitions in Function of g_k')
    xlabel('g_k')
    ylabel('NoReps')

    sgtitle('Steepest Descent Algorithm - Default g_k Value - Starting Point: ['+string(dg_sda{i}.StartingPoint(1))+','+string(dg_sda{i}.StartingPoint(2))+'] - \epsilon='+string(epsilon_sda)+' - '+f_label)
end

fprintf('Execution of Steepest Descent Algorithm - g_k="Minimize"\n')
ming_sda = {};

for i=1:size(starting_points,2)
    ming_sda{end+1} = struct('StartingPoint',starting_points(:,i));
    fprintf('Calculation of Minimum - StartingPoint=[%s]\n',strjoin(string(starting_points(:,i)),','))
    [ming_sda{i}.xValuesMin,ming_sda{i}.NoReps,ming_sda{i}.GamaList] = SteepestDescentAlgorithm(ming_sda{i}.StartingPoint,epsilon_sda,f,varsvect,"Minimize");
end

for i=2:length(ming_sda)
    figure()
    histogram(ming_sda{i}.GamaList)
    title('Steepest Descent Algorithm - g_k="Minimize" - Starting Point: ['+string(ming_sda{i}.StartingPoint(1))+','+string(ming_sda{i}.StartingPoint(2))+'] - \epsilon='+string(epsilon_sda)+' - '+f_label)
    subtitle('Variance of g_k Values Between Iterations - #g_k Values: '+string(length(ming_sda{i}.GamaList)))
    xlabel('g_k Intervals')
    ylabel('Amount of Values')
end

fprintf('Execution of Steepest Descent Algorithm - g_k="Armijo"\n')
armg_sda = {};
noarmreps = 15;
nobins = round(noarmreps/3);

for i=1:size(starting_points,2)
    armg_sda{end+1} = struct('StartingPoint',starting_points(:,i));
    mincalc = [];
    noreps = [];
    armparam = [];
    for j=1:noarmreps
        fprintf('Calculation of Minimum - StartingPoint=[%s] - #ArmijoIteration=%d\n',strjoin(string(starting_points(:,i)),','),j)
        [mincalc(end+1,:),noreps(end+1,:),~,armparam(end+1,:)] = SteepestDescentAlgorithm(armg_sda{i}.StartingPoint,epsilon_sda,f,varsvect,"Armijo","Random");
    end
    armg_sda{i}.xValuesMin = mincalc;
    armg_sda{i}.NoReps = noreps;
    armg_sda{i}.ArmijoParameters = armparam;
end


for i=2:length(armg_sda)
    figure()
    subplot(2,2,[1 3])
    hold on
    for j=1:size(armg_sda{i}.xValuesMin,1)
        plot(armg_sda{i}.xValuesMin(j,1),armg_sda{i}.xValuesMin(j,2),'o')
    end
    hold off
    title('Variance of (x_1,x_2) Values in Function of Armijo Parameters')
    xlabel('x_1 Values')
    ylabel('x_2 Values')
    labelnames = {};
    for j=1:size(armg_sda{i}.ArmijoParameters,1)
        labelnames{end+1} = strjoin(string(decimalAccuracyPoints(armg_sda{i}.ArmijoParameters(j,:),3)),',');
    end
    legend("alpha,beta,sigma = "+labelnames)
    
    f_temp_calc = [];
    for j=1:size(armg_sda{i}.xValuesMin,1)
        funmininparg = num2cell(transpose(armg_sda{i}.xValuesMin(j,:)));
        f_temp_calc(end+1) = feval(f_handle,funmininparg{:});
    end
    subplot(2,2,2)
    plot(1:1:noarmreps,f_temp_calc,'Marker','o','Color','Red')
    title('f(x_1,x_2) for Each Armijo Iteration')
    xlabel('NoIterations')
    ylabel('f(x_1,x_2)')

    subplot(2,2,4)
    plot(1:1:noarmreps,armg_sda{i}.NoReps,'Marker','o','Color','Blue')
    title('Algorithm''s Repetitions For Each Iteration')
    xlabel('No Armijo Iteration')
    ylabel('NoReps')

    sgtitle('Steepest Descent Algorithm - g_k="Armijo" & armijo_{parameters}="Random" - Starting Point: ['+string(armg_sda{i}.StartingPoint(1))+','+string(armg_sda{i}.StartingPoint(2))+'] - \epsilon='+string(epsilon_sda)+' - #ArmijoIterations: '+string(noarmreps)+' - '+f_label)

    figure()
    nexttile()
    histogram(armg_sda{i}.ArmijoParameters(:,1),nobins)
    title('Variance of "alpha" Armijo''s Parameter Between Iterations')
    xlabel('"alpha" Intervals')
    ylabel('Amount of Values')

    nexttile()
    histogram(armg_sda{i}.ArmijoParameters(:,2),nobins)
    title('Variance of "beta" Armijo''s Parameter Between Iterations')    
    xlabel('"beta" Intervals')
    ylabel('Amount of Values')

    nexttile()
    histogram(armg_sda{i}.ArmijoParameters(:,3),nobins)
    title('Variance of "sigma" Armijo''s Parameter Between Iterations')
    xlabel('"sigma" Intervals')
    ylabel('Amount of Values')

    sgtitle('Steepest Descent Algorithm - g_k="Armijo" & armijo_{parameters}="Random" - Starting Point: ['+string(armg_sda{i}.StartingPoint(1))+','+string(armg_sda{i}.StartingPoint(2))+'] - \epsilon='+string(epsilon_sda)+' - #ArmijoIterations: '+string(noarmreps)+' - '+f_label)
end