function [vmap] = sol(nx,ny,left_b, right_b, bot_b, top_b, d)
% -- Inputs
% nx = length of region
% ny = height of region
% left_b = left boundary
% right_b = right boundary
% bot_b = bottom boundary
% top_b = top boundary
% -- Ouputs
% vmap = potential map
g = sparse(nx*ny, nx*ny); % g matrix
f = zeros(nx*ny, 1); % f vector
if d == 'A'
    for i = 1:(nx)
        for j = 1:(ny)

            n = j + (i-1)*ny;

            if i == 1       
                g(n,:) = 0;
                g(n,n) = 1;
                f(n) = left_b;
            elseif i == nx
               g(n,:) = 0;
               g(n,n) = 1;
               f(n) = right_b;
             elseif j == 1
                nxm = j + (i - 2) * ny;
                nxp = j + i*ny;
                nyp = j + 1 + (i - 1) * ny;
    
                g(n, n) = -3;
                g(n, nxm) = 1;
                g(n, nxp) = 1;
                g(n, nyp) = 1;           
                 
             elseif j == ny
             
                nxm = j + (i-2)*ny;
                nxp = j + i*ny;
                nym = j-1 + (i-1)*ny;
    
                g(n,n) = -3;
                g(n, nxm) = 1;
                g(n, nxp) = 1;
                g(n, nym) = 1;
            else
                nxm = j + (i-2)*ny;
                nxp = j + i*ny;
                nym = j-1 + (i-1)*ny;
                nyp = j+1 + (i-1)*ny;

                g(n,n) = -4;

                g(n, nxm) = 1;
                g(n, nxp) = 1;
                g(n, nym) = 1;
                g(n, nyp) = 1;
            end
        end    
    end
elseif d == 'B'
    for i = 1:(nx)
        for j = 1:(ny)

            n = j + (i-1)*ny;

            if i == 1       
                g(n,:) = 0;
                g(n,n) = 1;
                f(n) = left_b; 
            elseif i == nx
               g(n,:) = 0;
               g(n,n) = 1;
               f(n) = right_b; 
             elseif j == 1
                g(n,:) = 0;
                g(n,n) = 1;      
                f(n) =  bot_b;                 
             elseif j == ny
                g(n,:) = 0;
                g(n,n) = 1;
                f(n) =  top_b;
             else
                nxm = j + (i-2)*ny;
                nxp = j + i*ny;
                nym = j-1 + (i-1)*ny;
                nyp = j+1 + (i-1)*ny;

                g(n,n) = -4;

                g(n, nxm) = 1;
                g(n, nxp) = 1;
                g(n, nym) = 1;
                g(n, nyp) = 1;
            end
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

% Part 1(b) Analytical series
if d == 'B'
    v_0 = 1;
    L = 3;
    W = 2;
    a = W;
    b = L;


    % Instead of going to infinity only go to a set number
    itr = 100; % 100 total iterations
    X = linspace(-b, b, nx);
    Y = linspace(0, a, ny);

    for r = 1:ny
        x(r,:) = X;
    end
    for c = 1:nx
        y(:,c) = Y;
    end
    
    soln = zeros(ny, nx);
    for i=1:itr
        n = 2*i - 1;
        soln = soln + (1./n)*((cosh((n.*pi.*x)./a))./(cosh((n.*pi.*b)./a))).*(sin((n.*pi.*y)./a));
    end

    series_soln = ((4.*v_0)./pi)*soln;

    figure(2);
    surf(series_soln);
    title("Series Solution of 1B")
end

end












