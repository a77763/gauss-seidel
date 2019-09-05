function [u,ITER]=PoissonGS(N, TOL)
% Finds the steady-state solution for the temperature distribution on a 
% square plate, for particular boundary conditions
% USES GAUSS-SEIDEL ITERATIONS
% N is the number of grid points in each direction (including boundary
% points)
% TOL is the tolerance for the stopping citeria
% IN THE GRID, USES i FOR ROW (ASCENDENT), j FOR COLUMN (LEFT-RIGHT)
% SEE Quinn's book, pp.330-332

% SET BOUNDARY VALUES 
% Temperature is zero at top and 100 on the other boundaries
u(1,1:N)=100;  % lower boundary 
u(1:N,1)=100;  % left boundary 
u(1:N,N)=100;  % right boundary 
u(N,1:N)=0;    % boundary above

w(1,1:N)=100;  % lower boundary 
w(1:N,1)=100;  % left boundary 
w(1:N,N)=100;  % right boundary 
w(N,1:N)=0;    % boundary above

% initial values for interior points
u(2:N-1,2:N-1)=50;

% COMPUTE STEADY STATE SOLUTION
DIFF=TOL+1;
ITER=0;
while DIFF>TOL
    for i=2:N-1
        for j=2:N-1
            w(i,j)=(w(i-1,j)+w(i,j-1)+u(i,j+1)+u(i+1,j))/4;
        end
    end
    DIFF=max(max(abs(w-u)));
    u=w; 
    ITER=ITER+1;
end

% DISPLAYS A COLOURFUL MAP OF DISTRIBUTIONS
u
u=flipud(u); image(u)






