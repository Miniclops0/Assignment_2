% Connor Warden
% 101078296

clc; close all; 

%These values can be easily changed to change to dimensions of every
%simulation. The set ratio is 3/2, so only nx needs to be changed.
nx = 75;        % Length
ny = nx*(2/3);  % Width

v_0 = 1;    % This is here incase it is needed for functions in the future,
            % the variable is named the same thing within the solutions

% Simulation type, A for 1A and B for 1B

d = 'A'; % if set to 'A' then 1A is performed

% Boundary conditions are set here, currently for Q 1A

left_b = v_0; % at x = 0, left
right_b = 0; % at x = L (max x), right. 
bot_b = 0; % at y = 0, bottom 
top_b = 0; % at y = max_y, top

[vmap] = sol(nx, ny, left_b, right_b, bot_b, top_b, d);
figure(3)
surf(vmap');

%To run both 1A and 1B at the same time the following two lines were added
right_b = v_0;
d = 'B';

[vmap] = sol(nx, ny, left_b, right_b, bot_b, top_b, d);
figure(5)
surf(vmap')

% Function call for question 2

sol_2(nx,ny);

