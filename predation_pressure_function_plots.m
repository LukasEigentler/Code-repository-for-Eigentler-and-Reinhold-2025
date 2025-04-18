%% predation pressure plot
% This script produces the appendix figure on the different predation
% pressures depending on the trait distribution of the prey.

clear;
close all;
alpha2 = 0.5;
c = linspace(-2,1/alpha2);
p = 0.5;
f2 = 1-alpha2*c;
mu_col = [0.5,1.5];
var = 0.2;
f = figure;
hold on

for bb = 1:length(mu_col)
    mu = mu_col(bb);
    x = 1/var*theta((c-mu)/var)/(phi((c(end)-mu)/var) - phi((c(1)-mu)/var));
    xp = trapz(c,f2.*x);
    pred_press = f2/(p+xp);
    graph(bb) = plot(c,pred_press);
%     plot(c,x)
end

% leg = legend(graph);
leg = legend(["Less well-defended prey", "More well-defended prey"]);
set(leg,"Interpreter", "latex")
xlabel("Defence trait, $c$", "Interpreter","latex")
ylabel("Per capita predation pressure, $\frac{f_2(c)}{p+X_p}$", "Interpreter","latex")
grid on

set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[10 5 10 11])

function theta = theta(x)
theta = 1/sqrt(2*pi)*exp(-x.^2/2);
end
function phi = phi(x)
phi = 1/2*(1+erf(x/sqrt(2)));
end