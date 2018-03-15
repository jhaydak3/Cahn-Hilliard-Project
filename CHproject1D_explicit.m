function CHproject1D_explicit 
close all
global D
global gamma
D = 100; 
gamma = .2;
dx = 5e-3;
dt = 1e-10;
save_fh = 'explicit_dt1e-10.mat';
[x1,~,c1] = my_CD(1e-10,1,dx,dt);
[x2,~,c2] = my_CD(1e-9,1,dx,dt);
[x3,~,c3] = my_CD(1e-8,1,dx,dt);
save(save_fh)
[x4,~,c4] = my_CD(1e-7,1,dx,dt);
save(save_fh)
[x5,~,c5] = my_CD(1e-6,1,dx,dt);
save(save_fh)
[x6,~,c6] = my_CD(1e-5,1,dx,dt);
save(save_fh)
[x7,~,c7] = my_CD(1e-4,1,dx,dt);
save(save_fh)
% hold all
% plot(x1,c1(:,end),'*-','MarkerSize',3)
% plot(x2,c2(:,end),'*-','MarkerSize',3)
% plot(x3,c3(:,end),'*-','MarkerSize',3)
%  plot(x4,c4(:,end),'*-','MarkerSize',3)
%   plot(x5,c5(:,end),'*-','MarkerSize',3)
% xlabel('t'); ylabel('c')
% title('Semi-implicit, D = 100, dx = 5e-3, dt = 1e-11, \gamma = .2')
end

function u = create_initial(u)
numx = length(u);
midpt = round(numx/3);
for i = 1:midpt
   %u(i) = sin(2*pi*(i/numx));
   u(i) = -1;
end
for i = midpt+1:midpt*2
  %u(i) = sin(2*pi*(i/numx));
  u(i) = 1;
end

for i = midpt*2+1:numx
   u(i) = -1; 
end

end
function [xout,tout,uout] = my_CD(t_f,x_f,dx,dt)
global D
global gamma
t = 0:dt:t_f;
x = 0:dx:x_f;
numx = length(x); numt = length(t);
uold  = zeros(numx,1);
uold = create_initial(uold);
unew = uold;
%plot initial conditions
% figure
plot(x,uold,'-*','MarkerSize',3)
xlabel('x'); ylabel('c');
title('Initial Conditions')

%mu = @(j,n) u(j,n)^3 - u(j,n) - gamma^2* (u(j+1,n)-2*u(j,n)+u(j-1,n))/dx^2;
for n = 1:numt-1
    display(['Iteration ' num2str(n) ' of '  num2str(numt)])
    for j = 3:numx-2
        %RHS = D* (mu(j+1,n) - 2*mu(j,n) - mu(j-1,n))/dx^2;
        %u(j,n+1) = RHS*dt + u(j,n);
        term1 = uold(j+1)^3 - 2*uold(j)^3 + uold(j-1)^3;
        term2 = uold(j+1) - 2*uold(j) + uold(j-1);
        term3= uold(j-2) - 4* uold(j-1) + 6*uold(j) - 4*uold(j+1) + uold(j+2);
        RHS = term1/dx^2 - term2/dx^2 - gamma^2*term3/dx^4;
        unew(j) = D*RHS*dt + uold(j);
    end
    j = 2;
    term1 = uold(j+1)^3 - 2*uold(j)^3 + uold(j-1)^3;
    term2 = uold(j+1) - 2*uold(j) + uold(j-1);
    term3= uold(1) - 4* uold(j-1) + 6*uold(j) - 4*uold(j+1) + uold(j+2);
    RHS = term1/dx^2 - term2/dx^2 - gamma^2*term3/dx^4;
    unew(j) = D*RHS*dt + uold(j);
    j = 1;
    unew(j) = unew(j+1);
    j = numx-1;
    term1 = uold(j+1)^3 - 2*uold(j)^3 + uold(j-1)^3;
    term2 = uold(j+1) - 2*uold(j) + uold(j-1);
    term3= uold(j-2) - 4* uold(j-1) + 6*uold(j) - 4*uold(j+1) + uold(numx);
    RHS = term1/dx^2 - term2/dx^2 - gamma^2*term3/dx^4;
    unew(j) = D*RHS*dt + uold(j);
    j = numx;
    unew(j) = unew(j-1);
    
    uold = unew;
    
end

xout = x;
tout = t;
uout = unew;
end