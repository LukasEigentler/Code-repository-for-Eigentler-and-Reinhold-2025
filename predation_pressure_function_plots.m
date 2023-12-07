%% predation pressure plot
% This script produces the appendix figure on the different predation
% pressures depending on the trait distribution of the prey.

clear;
close all;
cplot = 0.5;
c = linspace(0,1);
[~,cind] = min(abs(c-cplot));
p = 0.5;
alpha2 = 2;
f2 = alpha2-c;
b_col = [0.5,2];

f = figure;
hold on

for bb = 1:length(b_col)
    b = b_col(bb);
    x = 1-b/2 + b*c;
    xp = trapz(c,f2.*x);
    pred_press = f2/(p+xp);
    graph(bb) = plot(c,pred_press, "DisplayName","$F(X) = "+ num2str(1-b/2) + "+" + num2str(b)+"c$");
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