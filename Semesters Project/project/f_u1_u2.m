clc

u1_boundaries = [-1,2];
u2_boundaries = [-2,1];

u1 = u1_boundaries(1):0.01:u1_boundaries(2);
u2 = u2_boundaries(1):0.01:u2_boundaries(2);

[u1_grid,u2_grid] = meshgrid(u1,u2);

syms x1 x2
f_vars = [x1;x2];
f = sin(x1+x2)*sin(x2^2);
f_handle = matlabFunction(f,'Vars',f_vars);
f_val = f_handle(u1_grid,u2_grid);

figure()
surf(u1,u2,f_val)
title('Surface of f(u_1,u_2)=sin(u_1+u_2)*sin(u_2^2)')
xlabel('u_1')
ylabel('u_2')
zlabel('f(u_1,u_2)')
colorbar

figure()
contour(u1,u2,f_val,100)
title('Controur of f(u_1,u_2)=sin(u_1+u_2)*sin(u_2^2)')
subtitle('100 Isolines')
xlabel('u_1')
ylabel('u_2')
zlabel('f(u_1,u_2)')
colorbar