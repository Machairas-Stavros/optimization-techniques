%ConverganceDiagram - Some changes to the algorithm must be implemented in
%order to return the whole series of the (x_1,x_2) points.
a = [0.0342829557390931,0.221648315610874,0.417022211024702];
bbb = LevenbergMarquardtAlgorithm([-1;-1],10^-3,f,varsvect,"Armijo",a);
ccc = LevenbergMarquardtAlgorithm([1;1],10^-3,f,varsvect,0.95);
figure()
hold on
contour(x1_grid_values,x2_grid_values,y_values_f,100)
plot(bbb(1,:),bbb(2,:),'Color','Red','Marker','o')
plot(ccc(1,:),ccc(2,:),'Color','Blue','Marker','o')
hold off
title('Levenberg Marquardt Alogirthm - \epsilon=0.001')
subtitle('Method''s Convergence to Minimum Value')
ylabel('x_2 Values')
xlabel('x_1 Values')
legend('','g_k="Armijo" - [\alpha,\beta,s]=[0.034282955739093,0.2216483156108740.417022211024702] - Startin Point: [-1,-1]','Default g_k = 0.95 - Starting Point: [1,1]')