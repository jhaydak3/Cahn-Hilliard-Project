function CHproject_2D_implicit
close all
D = 100;
gamma = .2;
dr = .01;
dt = 1e-11;
t0 = 0;
tf1 = 1e-11;

save_str = '2d_implicit_rand_test2';
numx = length(0:dr:1);
c0 = create_initial_rand(numx,numx);
%imshow(c0)
tic
for i = 1:7
    istr = num2str(i);
    thisstr = ['[x' istr ',y' istr ',t' istr ',c' istr' ']'];
    if i == 1
        thisstr = [thisstr '=my_CD(0,tf1,1,1,dr,dr,dt,gamma,D,c0);'];
    else
        thisstr = [thisstr '=my_CD(t' num2str(i-1), ',t' num2str(i-1)];
        thisstr = [thisstr '*10,1,1,dr,dr,dt,gamma,D,c' num2str(i-1) ');'];
    end
    eval(thisstr)
    str = ['this_time' istr '=toc'];
    eval(str)
    save(save_str);
end

high = max(c(:));
low = min(c(:));
close all
imshow(c)
caxis([low high])
end
function A = create_initial_rand(numx,numy)
A = zeros(numx,numy);
for i = 1:numx
    for j = 1:numy
        A(i,j) = -1 + 2*rand;
    end
end
end
function A = create_initial_circ(numx,numy)
mdpntx = round(numx/2); mdpnty = round(numy/2);
A = zeros(numx,numy) -1;
for i = 1:numx
    for j = 1:numy
        dx = mdpntx - i;
        dy = mdpnty - j;
        if 2*dx^2 + dy^2 + 5*dx^2*dy - 20*dy^2*dx + 20*dx^3*dy - 10*dy^3*dx < 15
            A(i,j) = 1;
        end
        
    end
end
end

function [xout,yout,tout,new] = my_CD(t0,t_f,x_f,y_f,dx,dy,dt,gamma,D,u0)
t = t0:dt:t_f;
x = 0:dx:x_f;
y = 0:dy:y_f;
numx = length(x); numy = length(y); numt = length(t);
uold = u0;
unew = uold;
A = zeros(size(uold));
[nx,ny] = size(A);
b = zeros(nx,1);
rx = .5*D*dt/dx^2;
ry = .5*D*dt/dy^2;
rxx = .5*D*gamma^2*dt/dx^4;
ryy = .5*D*gamma^2*dt/dy^4;
for n = 1:numt-1
    %first do x direction
    for k = 1:numx %assume numx = numy
        for j = 1:numx
            if j == 2
                uoldjm2 = uold(k,1);
                uoldjm1 = uold(k,j-1);
            elseif j == 1
                uoldjm2 = uold(k,2);
                uoldjm1 = uold(k,1);
            else
                uoldjm1 = uold(k,j-1);
                uoldjm2 = uold(k,j-2);
            end
            if k == 2
                uoldkm2 = uold(1,j);
                uoldkm1 = uold(k-1,j);
            elseif k == 1
                uoldkm2 = uold(2,j);
                uoldkm1 = uold(1,j);
            else
                uoldkm1 = uold(k-1,j);
                uoldkm2 = uold(k-2,j);
            end
            if j == nx
                uoldjp1 = uold(k,nx);
                uoldjp2 = uold(k,nx-1);
            elseif j == nx-1
                uoldjp1 = uold(k,j+1);
                uoldjp2 = uold(k,nx);
            else
                uoldjp2 = uold(k,j+2);
                uoldjp1 = uold(k,j+1);
            end
            if k == nx
                uoldkp1 = uold(nx,j);
                uoldkp2 = uold(nx-1,j);
            elseif k == nx-1
                uoldkp1 = uold(k+1,j);
                uoldkp2 = uold(nx,j);
            else
                uoldkp1 = uold(k+1,j);
                uoldkp2 = uold(k+2,j);
            end
            if j > 1 && j < nx
                A(j,j-1) = rx - rx*(uoldjm1^2);
                A(j,j) = -2*rx + 2*rx*(uold(k,j)^2)+1;
                A(j,j+1) = rx - rx*(uoldjp1^2);
            elseif j==1
                A(j,j) = -rx + rx*uold(k,j)^2 + 1;
                A(j,j+1) = rx - rx*uoldjp1^2;
            elseif j == nx
                A(j,j) = -rx + rx*uold(k,j)^2 + 1;
                A(j,j-1) = rx - rx*uoldjm1^2;
            end
            b(j) = ry*(uoldkp1^3 - 2*uold(k,j)^3 + uoldkm1^3);
            b(j) = b(j) - ry*(uoldkp1 - 2*uold(k,j) + uoldkm1);
            b(j) = b(j) - rxx*(uoldjp2-4*uoldjp1+6*uold(k,j)-4*uoldjm1+uoldjm2);
            b(j) = b(j) - ryy*(uoldkp2-4*uoldkp1+6*uold(k,j)-4*uoldkm1+uoldkm2);
            b(j) = b(j) + uold(k,j);
        end
        this_sol = A \ b;
        unew(k,:) = this_sol;
    end
    uold = unew;
    A(:) = 0; b(:) = 0;%reset
    %now do y direction
    for j = 1:numx
        for k = 1:numx
            if j == 2
                uoldjm2 = uold(k,1);
                uoldjm1 = uold(k,j-1);
            elseif j == 1
                uoldjm2 = uold(k,2);
                uoldjm1 = uold(k,1);
            else
                uoldjm1 = uold(k,j-1);
                uoldjm2 = uold(k,j-2);
            end
            if k == 2
                uoldkm2 = uold(1,j);
                uoldkm1 = uold(k-1,j);
            elseif k == 1
                uoldkm2 = uold(2,j);
                uoldkm1 = uold(1,j);
            else
                uoldkm1 = uold(k-1,j);
                uoldkm2 = uold(k-2,j);
            end
            if j == nx
                uoldjp1 = uold(k,nx);
                uoldjp2 = uold(k,nx-1);
            elseif j == nx-1
                uoldjp1 = uold(k,j+1);
                uoldjp2 = uold(k,nx);
            else
                uoldjp2 = uold(k,j+2);
                uoldjp1 = uold(k,j+1);
            end
            if k == nx
                uoldkp1 = uold(nx,j);
                uoldkp2 = uold(nx-1,j);
            elseif k == nx-1
                uoldkp1 = uold(k+1,j);
                uoldkp2 = uold(nx,j);
            else
                uoldkp1 = uold(k+1,j);
                uoldkp2 = uold(k+2,j);
            end
            if k > 1 && k < nx
               A(k,k-1) = ry - ry*uoldkm1^2;
               A(k,k) = -2*ry + 2*ry*uold(k,j)^2 + 1;
               A(k,k+1) = ry - ry*uoldkp1^2;
            elseif k == 1
                A(k,k) = -ry + ry*uold(k,j)^2 + 1;
                A(k,k+1) = ry - ry*uoldkp1^2;
            elseif k ==nx 
                A(k,k) = -ry + ry*uold(k,j)^2 + 1;
                A(k,k-1) = ry - ry*uoldkm1^2;
            end
            b(k) = ry*(uoldjp1^3 - 2*uold(k,j)^3 + uoldjm1^3);
            b(k) = b(k) - rx*(uoldjp1 - 2*uold(k,j) + uoldjm1);
            b(k) = b(k) - rxx*( uoldjp2 -4*uoldjp1 + 6*uold(k,j) -4*uoldjm1 + uoldjm2 );
            b(k) = b(k) - ryy*(uoldkp2-4*uoldkp1 + 6*uold(k,j) -4*uoldkm1+ uoldkm2 );
            b(k) = b(k) + uold(k,j);
        end
        unew(:,j) = A \ b;
    end
    uold = unew;
end
xout = x;
yout = y;
tout = t(end);
new = unew;
end

