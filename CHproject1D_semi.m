function CHproject1D_semi
close all
global D
global gamma
D = 100; 
gamma = .2;
dx = 5e-3;
dt = 1e-10;
save_fh = 'semi_dt1e-10';
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
hold all
plot(x1,c1(:,end),'*-','MarkerSize',3)
plot(x2,c2(:,end),'*-','MarkerSize',3)
plot(x3,c3(:,end),'*-','MarkerSize',3)
plot(x4,c4(:,end),'*-','MarkerSize',3)
plot(x5,c5(:,end),'*-','MarkerSize',3)
plot(x6,c6(:,end),'*-','MarkerSize',3)
plot(x7,c7(:,end),'*-','MarkerSize',3)
xlabel('t'); ylabel('c')
legend('t=0','t=1e-10','t=1e-9','t=1e-8','t=1e-7','t=1e-6','t=1e-5','t=1e-4')
title('Semi-implicit, D = 100, dx = 5e-3, dt = 1e-11, \gamma = .2')
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

function [xout,tout,u] = my_CD(t_f,x_f,dx,dt)
global D
global gamma
t = 0:dt:t_f;
x = 0:dx:x_f;
numx = length(x); numt = length(t);
uold  = zeros(numx,1);
uold = create_initial(uold);
%plot initial conditions
figure
plot(x,uold,'-*','MarkerSize',3)
xlabel('x'); ylabel('c');
title('Initial Conditions')



%apply BC's
x = [x(1)-2*dx, x(1)-dx, x, x(end)+dx, x(end)+2*dx];
uold = [uold(2); uold(1); uold; uold(end); uold(end-1);];
new_numx = length(uold);
r = D*dt/dx^2;
for n = 1:numt-1
    display(['Iteration ' num2str(n) ' of '  num2str(numt)])
     A = zeros(new_numx,new_numx);
     b = zeros(new_numx,1);
     A(1,1) = 1;
     A(1,4) = -1;
     A(2,2) = 1;
     A(2,3) = -1;
     A(end,end) = 1; A(end,end-3) = -1;
     A(end-1,end-1) = 1; A(end-1,end-2) = -1;
    for j = 3:numx+2
        A(j,j-1) = r - r*(uold(j-1)^2);
        A(j,j) = -2*r + 2*r*(uold(j)^2)+1;
        A(j,j+1) = r - r*(uold(j+1)^2) ;
        b(j) = -gamma^2*D*dt/dx^4 * (uold(j+2)-4*uold(j+1) + 6*uold(j) -4*uold(j-1) +uold(j-2));
        b(j) = b(j) + uold(j);
    end
    unew = A \ b;
    uold = unew;
    
    
end

xout = x(3:numx+2);
tout = t;
u = unew(3:numx+2);
end