function [] = sol_2(nx,ny)
% -- Inputs
% nx = length of region
% ny = height of region

% Conductivity within the boxes
sigma_in = 0.01;

% cMap sets the conductivity within the boxes to sigma_in. The box size can
% be changed by altering the ratios within the if statements

cMap = ones(nx,ny);
for i = 1:nx
    for j = 1:ny
        if ((i>nx*0.45) && (i<nx*0.55) && (j<(ny*0.48)))
            cMap(i,j) = sigma_in;
        elseif ((i>nx*0.45) && (i<nx*0.55) && (j>(ny*0.52)))
            cMap(i,j) = sigma_in;
        end
    end
end

g = sparse(nx*ny); % g matrix
f = zeros(nx*ny, 1); % f vector

for i = 1:(nx)
    for j = 1:(ny)
        
        n = j + (i-1)*ny;
        
        if i == 1       
           g(n,n) = 1;
           f(n) = 1; 
        elseif i == nx
           g(n,n) = 1;
           f(n) = 0; 
        elseif j == 1
           nxm = j + (i - 2) * ny;
           nxp = j + (i) * ny;
           nyp = j + 1 + (i - 1) * ny;

           rxm = (cMap(i, j) + cMap(i - 1, j)) *0.5;
           rxp = (cMap(i, j) + cMap(i + 1, j)) *0.5;
           ryp = (cMap(i, j) + cMap(i, j + 1)) *0.5;

           g(n, n) = -(rxm+rxp+ryp);
           g(n, nxm) = rxm;
           g(n, nxp) = rxp;
           g(n, nyp) = ryp;
             
        elseif j == ny
           nxm = j + (i - 2) * ny;
           nxp = j + (i) * ny;
           nym = j - 1 + (i - 1) * ny;

           rxm = (cMap(i, j) + cMap(i - 1, j)) *0.5;
           rxp = (cMap(i, j) + cMap(i + 1, j)) *0.5;
           rym = (cMap(i, j) + cMap(i, j - 1)) *0.5;

           g(n, n) = -(rxm + rxp + rym);
           g(n, nxm) = rxm;
           g(n, nxp) = rxp;
           g(n, nym) = rym;
        else
           nxm = j + (i-2)*ny;
           nxp = j + (i)*ny;
           nym = j-1 + (i-1)*ny;
           nyp = j+1 + (i-1)*ny;

           rxm = (cMap(i,j) + cMap(i-1,j))*0.5;
           rxp = (cMap(i,j) + cMap(i+1,j))*0.5;
           rym = (cMap(i,j) + cMap(i,j-1))*0.5;
           ryp = (cMap(i,j) + cMap(i,j+1))*0.5;

           g(n,n) = -(rxm+rxp+rym+ryp);
           g(n,nxm) = rxm;
           g(n,nxp) = rxp;
           g(n,nym) = rym;
           g(n,nyp) = ryp;
        end
    end    
end

v = g\f;


vmap = zeros(nx,ny);
for i = 1:(nx)
    for j = 1:(ny)
        n = j + (i-1) * ny;
        vmap(i, j) = v(n);
    end
end

for i = 1:(nx)
    for j = 1:(ny)
        if i == 1
            Ex(i,j) = vmap(i+1,j) - vmap(i,j);
        elseif i == nx
            Ex(i,j) = vmap(i,j) - vmap(i-1,j);
        else
            Ex(i,j) = (vmap(i+1,j) - vmap(i-1,j))*(0.5);
        end
        if j == 1
            Ey(i,j) = vmap(i,j+1) - vmap(i,j);
        elseif j == ny
            Ey(i,j) = vmap(i,j) - vmap(i,j-1);
        else
            Ey(i,j) = (vmap(i,j+1) - vmap(i,j-1))*(0.5);
        end
    end
end

Ex = -Ex;
Ey = -Ey;

i_x = cMap .* Ex;
i_y = cMap .* Ey;


figure(4)
surf(vmap')
title("Potential of 2")

figure(5)
surf(cMap')
title("Conductivity Plot")

figure(6)
quiver(Ex', Ey');
title("Electric Field")

figure(7)
quiver(i_x', i_y')
title("Current Density")

end



