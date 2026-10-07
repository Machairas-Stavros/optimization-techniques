clear
clc

rng('default')

u1_inter = [-1,2];
u2_inter = [-2,1];

u1 = u1_inter(1):0.01:u1_inter(2);
u2 = u2_inter(1):0.01:u2_inter(2);

[u1_grid,u2_grid] = meshgrid(u1,u2);

m1_inter = [0,1];
m2_inter = [-1,0];
s1_inter = [0.01,1.5];
s2_inter = [0.01,1.5];
c_inter = [-1,1];

nb = 5;
cb = 3;
ng = 15;
nsp = 500;

m1 = m1_inter(1):(m1_inter(2)-m1_inter(1))/(2^nb-1):m1_inter(2);
m2 = m2_inter(1):(m2_inter(2)-m2_inter(1))/(2^nb-1):m2_inter(2);
s1 = s1_inter(1):(s1_inter(2)-s1_inter(1))/(2^nb-1):s1_inter(2);
s2 = s2_inter(1):(s2_inter(2)-s2_inter(1))/(2^nb-1):s2_inter(2);
c = c_inter(1):(c_inter(2)-c_inter(1))/(2^nb-1):c_inter(2);

epsilon = 0.08;
dap = 5;
max_reps = 50;
total_iter = 100;

initial_population = {};

for i=1:total_iter
    pop_mse = [];
    min_pop_mse = inf;
    iter = 0;
    population = createPopulation(m1,m2,s1,s2,c,ng,nsp);
    chromosome = [];
    while true
        iter = iter+1;
        [population,pop_mse(:,end+1),enditer,iter_chromosome,min_pop_mse,stpr] = generationUpdateCycle(population,m1',m2',s1',s2',c',u1_grid,u2_grid,max_reps,iter,epsilon,min_pop_mse,nb,cb,dap);
        fprintf('%d Split - %d Iteration - MSE Min Value: %.3f - MSE Mean Value: %.3f - MSE Max Value: %.3f\n',i,iter,min(pop_mse(:,iter)),mean(pop_mse(:,iter),"omitmissing"),max(pop_mse(:,iter)));

        if ~isempty(iter_chromosome) && (isempty(stpr) || stpr=="E")
            chromosome = iter_chromosome;
        end

        if enditer
            break
        end
    end
    initial_population{end+1} = struct('Chromosome',chromosome,'MSE',pop_mse,'MinMSE',min_pop_mse);
end

clear pop_mse enditer chromosome iter iter_chromosome min_pop_mse stpr population

ultimate_population = [];

for i=1:length(initial_population)
    ultimate_population(:,:,end+1) = initial_population{i}.Chromosome;
end

ultimate_population = ultimate_population(:,:,randsample(size(ultimate_population,3),500,true));

pop_mse = [];
iter = 0;
min_pop_mse = inf;
ultimate_chromosome = [];

max_reps = 500;
epsilon = 0.01;

while true
    iter = iter+1;
    [ultimate_population,pop_mse(:,end+1),enditer,iter_chromosome,min_pop_mse,stpr] = generationUpdateCycle(ultimate_population,m1',m2',s1',s2',c',u1_grid,u2_grid,max_reps,iter,epsilon,min_pop_mse,nb,cb,dap);
    fprintf('%d Iteration - MSE Min Value: %.3f - MSE Mean Value: %.3f - MSE Max Value: %.3f\n',iter,min(pop_mse(:,iter)),mean(pop_mse(:,iter),"omitmissing"),max(pop_mse(:,iter)));

    if ~isempty(iter_chromosome) && (isempty(stpr) || stpr=="E")
        ultimate_chromosome = iter_chromosome;
    end

    if enditer
        break
    end
end
ultimate_population = struct('Chromosome',ultimate_chromosome,'MSE',pop_mse,'MinMSE',min_pop_mse);

clear pop_mse enditer iter ultimate_chromosome iter_chromosome min_pop_mse stpr

clf
close all

initial_population_mse = [];

for i=1:length(initial_population)
    initial_population_mse(end+1) = initial_population{i}.MinMSE;
end

[min_mse_ip_val,min_mse_ip_ind] = min(initial_population_mse);

figure()
subplot(2,1,1)
histogram(initial_population_mse,25)
title('Histogram of Minimum MSE Achieved from 100 Different Initial Populations')
subtitle('15 Gaussians Used')
xlabel('Minimum MSE Value')
ylabel('Times Appeared')

subplot(2,1,2)
hold on
plot(initial_population_mse,'Marker','o')
plot(min_mse_ip_ind,min_mse_ip_val,'Marker','o','MarkerFaceColor','red')
hold off
title('Minimum MSE Achieved Over Each Different Initial Population')
subtitle('15 Gaussians Used')
xlabel('No Iteration')
ylabel('Minimum MSE Value')
legend('Time Series','Overall Minimum Value')

clear initial_population_mse min_mse_ip_val min_mse_ip_ind

min_mse_up = min(ultimate_population.MSE);
[min_mse_up_val,min_mse_up_ind] = min(min_mse_up);

figure()
subplot(2,1,1)
histogram(min(ultimate_population.MSE),25)
title('Histogram of Minimum MSE Achieved Over Each Iteration from Ultimate Population')
subtitle('15 Gaussians Used')
xlabel('Minimum MSE Value')
ylabel('Times Appeared')

subplot(2,1,2)
hold on
plot(min_mse_up,'Marker','o')
plot(min_mse_up_ind,min_mse_up_val,'Marker','o','MarkerFaceColor','red')
hold off
title('Minimum MSE Achieved Over Each Iteration from Ultimate Population')
subtitle('15 Gaussians Used')
xlabel('No Iteration')
ylabel('Minimum MSE Value')
legend('Time Series','Overall Minimum Value')

clear min_mse_up min_mse_up_val min_mse_up_ind

initial_population_mse = [];

for i=1:length(initial_population)
    initial_population_mse(end+1) = initial_population{i}.MinMSE;
end

[~,initial_population_mse_ind] = min(initial_population_mse);

inpsurf = calculateGaussianLinearCombo(u1_grid,u2_grid,initial_population{initial_population_mse_ind}.Chromosome);
ulpsurf = calculateGaussianLinearCombo(u1_grid,u2_grid,ultimate_population.Chromosome);
evalsurf = evaluationFunct(u1_grid,u2_grid);

figure()
surf(u1_grid,u2_grid,inpsurf)
title('Estimated Surface of the Best Chromosome from Initial Population')
subtitle('15 Gaussians Used')
xlabel('u_1')
ylabel('u_2')
zlabel('Estimated Surface')
colorbar

figure()
contour(u1_grid,u2_grid,sqrt((inpsurf-evalsurf).^2),100)
title('Contour of Squared Error at Each Point - Best Chromosome from Initial Population - 100 Isolines')
subtitle('15 Gaussians Used')
xlabel('u_1')
ylabel('u_2')
colorbar

figure()
surf(u1_grid,u2_grid,ulpsurf)
title('Estimated Surface of the Chromosome from Ultimate Population')
subtitle('15 Gaussians Used')
xlabel('u_1')
ylabel('u_2')
zlabel('Estimated Surface')
colorbar

figure()
contour(u1_grid,u2_grid,sqrt((ulpsurf-evalsurf).^2),100)
title('Contour of Squared Error at Each Point - Chromosome from Ultimate Population - 100 Isolines')
subtitle('15 Gaussians Used')
xlabel('u_1')
ylabel('u_2')
colorbar

clear initial_population_mse initial_population_mse_ind inpsurf ulpsurf evalsurf

u1_test = u1_inter(1)+0.005:0.01:u1_inter(2)-0.005;
u2_test = u2_inter(1)+0.005:0.01:u2_inter(2)-0.005;

[u1_grid_test,u2_grid_test] = meshgrid(u1_test,u2_test);

initial_population_mse = [];
for i=1:length(initial_population)
    initial_population_mse(end+1) = initial_population{i}.MinMSE;
end

[~,initial_population_mse_ind] = min(initial_population_mse);

inpsurf_test = calculateGaussianLinearCombo(u1_grid_test,u2_grid_test,initial_population{initial_population_mse_ind}.Chromosome);
ulpsurf_test = calculateGaussianLinearCombo(u1_grid_test,u2_grid_test,ultimate_population.Chromosome);
evalsurf_test = evaluationFunct(u1_grid_test,u2_grid_test);

mse_ip_test = mse3D(inpsurf_test,evalsurf_test);
mse_up_test = mse3D(ulpsurf_test,evalsurf_test);

clear u1_test u2_test u1_grid_test u2_grid_test initial_population_mse initial_population_mse_ind inpsurf_test ulpsurf_test evalsurf_test