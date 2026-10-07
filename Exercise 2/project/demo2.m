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

starting_point = [-1;-0.8];
gama_test_values = 0.05:0.05:1;

epsilon_na = 10^-10;

fprintf('Execution of Newton Algorithm - Default g_k\n')

dg_na = struct('StartingPoint',starting_point);
mincalc = [];
noreps = [];
for i=1:length(gama_test_values)
    fprintf('Calculation of Minimum - StartingPoint=[%s] - g_k=%.2f\n',strjoin(string(starting_point),','),gama_test_values(i))
    [mincalc(end+1,:),noreps(end+1,:)] = NewtonAlgorithm(dg_na.StartingPoint,epsilon_na,f,varsvect,gama_test_values(i));
end
dg_na.xValuesMin = mincalc;
dg_na.NoReps = noreps;

figure()
subplot(2,2,[1 3])
hold on
for i=1:size(dg_na.xValuesMin,1)
    plot(dg_na.xValuesMin(i,1),dg_na.xValuesMin(i,2),'o')
end
hold off
title('Variance of (x_1,x_2) Values in Function of g_k')
xlabel('x_1 Values')
ylabel('x_2 Values')
legend("gk = "+num2cell(string(gama_test_values)))

f_temp_calc = [];
for i=1:size(dg_na.xValuesMin,1)
    funmininparg = num2cell(transpose(dg_na.xValuesMin(i,:)));
    f_temp_calc(end+1) = feval(f_handle,funmininparg{:});
end

subplot(2,2,2)
plot(gama_test_values,f_temp_calc,'Marker','o','Color','Red')
title('f(x_1,x_2) in Function of g_k')
xlabel('g_k')
ylabel('f(x_1,x_2)')

subplot(2,2,4)
plot(gama_test_values,dg_na.NoReps,'Marker','o','Color','Blue')
title('Algorithm''s Repetitions in Function of g_k')
xlabel('g_k')
ylabel('NoReps')

sgtitle('Newton Algorithm - Default g_k Value - Starting Point: ['+string(dg_na.StartingPoint(1))+','+string(dg_na.StartingPoint(2))+'] - \epsilon='+string(epsilon_na)+' - '+f_label)

fprintf('Execution of Newton Algorithm - g_k="Minimize"\n')

ming_na = struct('StartingPoint',starting_point);
fprintf('Calculation of Minimum - StartingPoint=[%s]\n',strjoin(string(starting_point),','))
[ming_na.xValuesMin,ming_na.NoReps,ming_na.GamaList] = NewtonAlgorithm(ming_na.StartingPoint,epsilon_na,f,varsvect,"Minimize");

figure()
histogram(ming_na.GamaList)
title('Newton Algorithm - g_k="Minimize" - Starting Point: ['+string(ming_na.StartingPoint(1))+','+string(ming_na.StartingPoint(2))+'] - \epsilon='+string(epsilon_na)+' - '+f_label)
subtitle('Variance of g_k Values Between Iterations - #g_k Values: '+string(length(ming_na.GamaList)))
xlabel('g_k Intervals')
ylabel('Amount of Values')

fprintf('Execution of Newton Algorithm - g_k="Armijo"\n')
noarmreps = 15;
nobins = round(noarmreps/3);

armg_na = struct('StartingPoint',starting_point);
mincalc = [];
noreps = [];
armparam = [];
for i=1:noarmreps
    fprintf('Calculation of Minimum - StartingPoint=[%s] - #ArmijoIteration=%d\n',strjoin(string(starting_point),','),i)
    [mincalc(end+1,:),noreps(end+1,:),~,armparam(end+1,:)] = NewtonAlgorithm(armg_na.StartingPoint,epsilon_na,f,varsvect,"Armijo","Random");
end
armg_na.xValuesMin = mincalc;
armg_na.NoReps = noreps;
armg_na.ArmijoParameters = armparam;

figure()
subplot(2,2,[1 3])
hold on
for i=1:size(armg_na.xValuesMin,1)
    plot(armg_na.xValuesMin(i,1),armg_na.xValuesMin(i,2),'o')
end
hold off
title('Variance of (x_1,x_2) Values in Function of Armijo Parameters')
xlabel('x_1 Values')
ylabel('x_2 Values')
labelnames = {};
for i=1:size(armg_na.ArmijoParameters,1)
    labelnames{end+1} = strjoin(string(decimalAccuracyPoints(armg_na.ArmijoParameters(i,:),3)),',');
end
legend("alpha,beta,sigma = "+labelnames)
    
f_temp_calc = [];
for i=1:size(armg_na.xValuesMin,1)
    funmininparg = num2cell(transpose(armg_na.xValuesMin(i,:)));
    f_temp_calc(end+1) = feval(f_handle,funmininparg{:});
end
subplot(2,2,2)
plot(1:1:noarmreps,f_temp_calc,'Marker','o','Color','Red')
title('f(x_1,x_2) for Each Armijo Iteration')
xlabel('NoIterations')
ylabel('f(x_1,x_2)')

subplot(2,2,4)
plot(1:1:noarmreps,armg_na.NoReps,'Marker','o','Color','Blue')
title('Algorithm''s Repetitions For Each Iteration')
xlabel('No Armijo Iteration')
ylabel('NoReps')

sgtitle('Nweton Algorithm - g_k="Armijo" & armijo_{parameters}="Random" - Starting Point: ['+string(armg_na.StartingPoint(1))+','+string(armg_na.StartingPoint(2))+'] - \epsilon='+string(epsilon_na)+' - #ArmijoIterations: '+string(noarmreps)+' - '+f_label)

figure()
nexttile()
histogram(armg_na.ArmijoParameters(:,1),nobins)
title('Variance of "alpha" Armijo''s Parameter Between Iterations')
xlabel('"alpha" Intervals')
ylabel('Amount of Values')

nexttile()
histogram(armg_na.ArmijoParameters(:,2),nobins)
title('Variance of "beta" Armijo''s Parameter Between Iterations')    
xlabel('"beta" Intervals')
ylabel('Amount of Values')

nexttile()
histogram(armg_na.ArmijoParameters(:,3),nobins)
title('Variance of "sigma" Armijo''s Parameter Between Iterations')
xlabel('"sigma" Intervals')
ylabel('Amount of Values')

sgtitle('Nweton Algorithm - g_k="Armijo" & armijo_{parameters}="Random" - Starting Point: ['+string(armg_na.StartingPoint(1))+','+string(armg_na.StartingPoint(2))+'] - \epsilon='+string(epsilon_na)+' - #ArmijoIterations: '+string(noarmreps)+' - '+f_label)