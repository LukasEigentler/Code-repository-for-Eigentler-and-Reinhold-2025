function dvdt = pred_prey_prey_defence_ode(v,c,M,d,alpha1,alpha2,ph,gamma,m2,m1,alt,s)
% This function contains the ODEs fed into the Matlab ODE solver
dvdt = zeros(M+1,1); % M entries for X, 1 entry for total Y
X =  v(1:M)';
Y = v(end)';

if M>1
    if alt < 3
        relpred = max(zeros(1,length(c)),(1-alpha2*c)).*X;
    elseif alt == 3
        relpred = max(zeros(1,length(c)),(1-(1-s)*c.^2./(c.^2+alpha2^(-2)))).*X;
    else
        error("Wrong value for alt")
    end
    xbar = trapz(c,relpred);
    dc = c(2)-c(1);
    Xtot = sum(X)*dc;
else
    xbar = (2-alpha2)*X\2;
    dc = NaN;
    Xtot = X;
end
if alt < 3
    p = max(zeros(1,length(c)),(1-alpha2*c))*Y./(ph+xbar);
elseif alt == 3
    p = max(zeros(1,length(c)),(1-(1-s)*c.^2./(c.^2+alpha2^(-2))))*Y./(ph+xbar);
end
if alt == 0
    g = (1-alpha1*c).*(1-Xtot)-m1;
elseif alt == 1 || alt == 3
    g = (1-alpha1*c)-Xtot;
elseif alt == 2
    g = (1-alpha1*c.^2)-Xtot;
else
    error("Wrong value for alt")
end



ccu = 2:M-1;

if M>1
    dvdt(ccu) = X(ccu).*(g(ccu) - p(ccu)) + d*(X(ccu+1)-2*X(ccu)+X(ccu-1))/(dc^2); % prey
    dvdt(1) = X(1)*(g(1) - p(1)) + d*(2*X(2)-2*X(1))/(dc^2);
    dvdt(M) = X(M)*(g(M) - p(M)) + d*(2*X(M-1) - 2*X(M))/(dc^2);
else
    dvdt(1) = X(1)*(g(1) - p(1));
end
dvdt(M+1) = gamma*xbar*Y/(ph+xbar) - m2*Y; % predator


