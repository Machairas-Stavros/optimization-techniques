clf
close all
clear
clc

syms x1 x2
varsvect = [x1;x2];

f = x1^2/3+3*x2^2;
f_handle = matlabFunction(f,'Vars',{'x1','x2'});
f_label = "f(x_1,x_2)=x_1^2/3+3*x_2^2";
x_values = -10:0.1:10;

[x1_grid_values,x2_grid_values] = meshgrid(x_values,x_values);
y_values_f = f_handle(x1_grid_values,x2_grid_values);

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

starting_points = [5,-5,8,;-5,10,-10];
epsilon = 0.01;
sk = [5,15,0.1];
gk = [0.5,0.1,0.2];

restriction_boundaries = [-10,5;-8,12];

dg_sda = {};
for i=1:size(starting_points,2)
    dg_sda{end+1} = struct('StartingPoint',starting_points(:,i),'Gama',gk(i),'Sk',sk(i));
    [dg_sda{i}.xValuesMin,dg_sda{i}.yValuesMin,dg_sda{i}.NoReps] = SteepestDescentAlgorithmProjection(dg_sda{i}.StartingPoint,epsilon,dg_sda{i}.Gama,dg_sda{i}.Sk,f,varsvect,restriction_boundaries);
end

for i=1:size(starting_points,2)
    figure()
    nexttile()
    hold on
    contour(x1_grid_values,x2_grid_values,y_values_f,100)
    plot(dg_sda{i}.xValuesMin(1,:),dg_sda{i}.xValuesMin(2,:),'Marker','o','Color','Red')
    hold off
    title('Convergence of (x1,x2)')
    subtitle('Contour - 100 Isolines')
    xlabel('x_1 Values')
    ylabel('x_2 Values')
    legend(horzcat("contour","gk = "+string(gk(i))))
    colorbar
    
    nexttile()
    plot(dg_sda{i}.yValuesMin,'Marker','o','Color','Blue')
    title('Converge of f(x_1,x_2)')
    xlabel('No. Iterations')
    ylabel('f(x_1,x_2)')
    
    sgtitle('Steepest Descent Projection Algorithm - g_k='+string(gk(i))+' - s_k='+string(sk(i))+' - Starting Point: ['+string(dg_sda{i}.StartingPoint(1))+','+string(dg_sda{i}.StartingPoint(2))+'] - \epsilon='+string(epsilon)+' - '+f_label)
end