%% Predator-prey model with prey defence trait: single run
% This script solves the model once and plots the output. Also plots Fig 5
% of paper, which compares the fitness of single genotypes across 
% different simulations if corresponding parameter values are chosen.

clear; 
close all;

%% Parameters
d=1e-3; % mutation rate
alpha1 = 0.75; % cost of prey defence
alpha2 = 0.5; % prey defence efficiency
m1 = 0.2; %prey mortality
m2 = 0.2; %pred mortality (LV only)
ph = 0.5; %predation half saturation constant (extension only)
gamma = 4; % prey to predator conversion

%% Mesh
cmax = 1; %Space domain size of interest
tmax = 1000; %Integration range for solver
M = 2^8; %Number of trait points
% option to reduce model to a non-evolutionary model
if M>1
    c=linspace(0,cmax,M);
else
    c=0.0;
end

%% initial condition
u0 = 0.5*ones(1,length(c)); % prey 
u0(end+1) = 0.5; % predator


%% ODE Solver
options = odeset('Stats', 'off','MaxStep',1e-2,'NonNegative',1:M+1); 
[t,v,totalprey,medianc,meanc,L,v_op,totalprey_op,t_op,medianc_op,meanc_op,varc_op,phaselag_prey_pred,phaselag_pred_trait,phaselag_mean_var,varc] = prey_defence_single_run_fun(c,M,d,alpha1,alpha2,ph,gamma,m2,m1,tmax,u0,options);

%% visualisations 

plotstart = t(end) - 2*t_op(end);
% plotstart = 00;
figure;
subplot(2,2,1)
plot(t,totalprey)
xlabel("t")
ylabel("Prey")
% xlim([plotstart,tmax])

subplot(2,2,2)
plot(t,v(:,M+1))
xlabel("t")
ylabel("Pred")
% xlim([plotstart,tmax])

subplot(2,2,3)
plot(t,totalprey)
hold on
plot(t,v(:,M+1))
xlabel("Time, t")
ylabel("both")
yyaxis right
plot(t,medianc, '--')
hold on
plot(t,meanc, '--')
ylabel("Mean trait")
% xlim([plotstart,tmax])

subplot(2,2,4)
T = length(t);
plot(totalprey(1:ceil(9/10*T)),v(1:ceil(9/10*T),M+1))
hold on
plot(totalprey(ceil(9/10*T):end),v(ceil(9/10*T):end,M+1))
xlabel("Prey")
ylabel("Pred")

%% plot prey pred densities and mean trait
f = figure;
plot(t_op,totalprey_op)
hold on
% grid on
plot(t_op,v_op(:,M+1))
% xlabel("Time, t")
xticklabels([])
ylabel("Prey and pred. dens.")
legend("Total prey", "Predator")
grid on
% ylim([0,1.1])
% yyaxis right
% plot(t_op,meanc_op, '--')
% ylabel("Mean trait")
% % xlim([plotstart,t(end)])
% ylim([0,1])

set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',10)
set(f,'Units','centimeters')
set(f,'Position',[18 1 17/2 4])
% exportgraphics(f,"../../Ecol_paper/figures/example_po_2.eps","Resolution",2000)
%% calc IQR
% np = length(c);
%     interq_trait=zeros(1,length(t)); traitdist = zeros(np,length(t));
%     for i=1:length(t) % loop through all times
%         traitdist(1:np,i) = v(i,1:np)/sum(v(i,1:np)); % frequency of trait value
%         csum = cumsum(traitdist(:,i)); %cdf of trait dist
%         cspace = linspace(0,1,1000); % larger c vector
%         csum = interp1(c',csum,cspace); % interpolate onto larger c vector
%         p25ind = find(csum<0.25); % find 25 percentile
%         if ~isempty(p25ind) 
%             p25ind = p25ind(end)+1; cp25 = cspace(p25ind);
%         else
%             cp25 = cspace(1);
%         end
%         p75ind = find(csum<0.75); % find 75 percentile
%         if ~isempty(p75ind)
%             p75ind = p75ind(end); cp75 = cspace(p75ind);
%         else
%             cp75 = cspace(end);
%         end
%         interq_trait(i) = sqrt(varc(i)); % interquartile range of traits
%     end
   
std_trait = sqrt(varc);
%% plot mean and IQR
f2 = figure;
col = lines;
hold on
grid on
xlabel("Time, t")
plot(t_op,meanc_op, 'color', 'k')
plot(t_op,sqrt(varc_op),'-','color',col(4,:))
try
plot([t_op(islocalmax(meanc_op)),t_op(islocalmax(meanc_op))],[0,meanc_op(islocalmax(meanc_op))],'-','color', [.5 .5 .5],'LineWidth',0.5)
plot([t_op(islocalmin(meanc_op)),t_op(islocalmin(meanc_op))],[0,meanc_op(islocalmin(meanc_op))],'-','color', [.5 .5 .5],'LineWidth',0.5)
end
ylabel("Trait mean and std")
ylim([0,1])
leg = legend("Mean", "Std");
leg.FontSize = 9;
% xlim([plotstart,t(end)])

set(f2,'Windowstyle','normal')
set(findall(f2,'-property','FontSize'),'FontSize',10)
set(f2,'Units','centimeters')
set(f2,'Position',[18 1 8.5 4])

% single time points over space

figure
subplot(2,2,1)
plot(c,v(1,1:M))
xlabel("c")
ylabel("Prey")

subplot(2,2,2)
plot(c,v(1,M+1)*ones(1,length(c)))
xlabel("c")
ylabel("Pred")

% temp_ind = find(islocalmax(totalprey));
% plotind_single = temp_ind(end-1);

temp_ind = find(islocalmax(v(:,M+1)));
plotind_single = length(v(:,1)); %temp_ind(end-1);

subplot(2,2,3)
plot(c,v(plotind_single,1:M))
xlabel("c")
ylabel("Prey")

subplot(2,2,4)
plot(c,v(plotind_single,M+1)*ones(1,length(c)))
xlabel("c")
ylabel("Pred")

%% distribution plot single time point
try
f3 = figure;
% plotind_single = length(t); % end point of simulation
% tempind = find(islocalmin(meanc_op)); % at min meanc
tempind = find(islocalmax(meanc_op)); % at max meanc
% tempind = find(t_op<7);
plotind_single =  tempind(end);
plot(c,v_op(plotind_single,1:M)/((c(2)-c(1))*sum(v_op(plotind_single,1:M))))
% xlabel("Trait value, c")
% ylabel("Prey frequency")
grid on
ylim([0,4.5])

set(f3,'Windowstyle','normal')
set(findall(f3,'-property','FontSize'),'FontSize',11)
set(f3,'Units','centimeters')
set(f3,'Position',[18 1 17/4 17/4])


% saveas(f,"Paper/figures/example_po_f1_eq_f2", 'epsc')

end

%%

if max(abs(max(totalprey_op)-min(totalprey_op))./mean(totalprey_op)) > 1e-2 && max(totalprey_op) > 1e-2
    disp("Limit cycle solution detected /n")
    const = 0;
    % xlim([plotstart,tmax])
    % save("start_data_orbit", "v_op", "c", "t_op")
else
    disp("Constant solution detected /n")
    v_op = max(v_op);
    const = 1;
    
    
end
%%
f1 = figure;
subplot(2,3,1)
plot(t_op,totalprey_op)
hold on
plot(t_op,v_op(:,M+1))
xlabel("t")
ylabel("both")
yyaxis right
plot(t_op,meanc_op, '--')
plot(t_op,varc_op,'--','color','k')
ylabel("Mean trait/std")

subplot(2,3,2)
[maxtrait,maxtraitind] = max(meanc_op);
if const == 1
    plot(c,v_op(1,1:M))
else
    plot(c,v_op(maxtraitind,1:M))
end
xlabel("c")
ylabel("Prey")
grid on
title("Trait distribution at max median trait")

try
subplot(2,3,3)
[mintrait,mintraitind] = min(meanc_op);
plot(c,v_op(mintraitind,1:M))
xlabel("c")
ylabel("Prey")
grid on
title("Trait distribution at min median trait")
end

try
subplot(2,3,4)
minmaxhalftrait = mintrait + 0.5*(maxtrait - mintrait); % diff between max and min median trait
tempind = find(meanc_op>minmaxhalftrait);
plot(c,v_op(tempind(1),1:M))
xlabel("c")
ylabel("Prey")
grid on
title("Trait distribution at midpoint (inreasing)")
end
try
subplot(2,3,5)
plot(c,v_op(tempind(end),1:M))
xlabel("c")
ylabel("Prey")
grid on
title("Trait distribution at midpoint (decreasing)")
end
fprintf("Solution properties: "...
    + "\n Prey biomass (mean): "+ num2str(mean(totalprey_op),'%4.2f') + ", Prey biomass (max): "+ num2str(max(totalprey_op),'%4.2f') + ", Prey biomass (min): "+ num2str(min(totalprey_op),'%4.2f') ...
    + "\n Pred biomass (mean): "+ num2str(mean(v_op(:,end)),'%4.2f') + ", Prey biomass (max): "+ num2str(max(v_op(:,end)),'%4.2f') + ", Prey biomass (min): "+ num2str(min(v_op(:,end)),'%4.2f') ... 
    + "\n Prey mean trait (mean): "+ num2str(mean(meanc_op),'%4.2f') + ", Prey mean trait (max): "+ num2str(max(meanc_op),'%4.2f') + ", Prey mean trait (min): "+ num2str(min(meanc_op),'%4.2f') ...
    + "\n Prey trait iqr (mean): "+ num2str(mean(varc_op),'%4.2f') + ", Prey trait iqr (max): "+ num2str(max(varc_op),'%4.2f') + ", Prey trait iqr (min): "+ num2str(min(varc_op),'%4.2f') ...
    + "\n Wavelength: " + num2str(L,'%4.2f') ...
    + "\n phaseshift prey-pred: " + num2str(phaselag_prey_pred,'%4.2f') + ", phaseshift pred-trait: " + num2str(phaselag_pred_trait,'%4.2f')+"\n")


%% calc prey fitness for each genotype at one time point
if ~(max(abs(max(totalprey_op)-min(totalprey_op))./mean(totalprey_op)) > 1e-2 && max(totalprey_op) > 1e-2) % only for constant sol
    X = v_op(1:end-1); Y = v_op(end);
else
    [~,maxiqrind] = max(varc_op);
    X = v_op(maxiqrind,1:end-1); Y = v_op(maxiqrind,end);
end
relpred = (1-alpha2*c).*X;
xbar = trapz(c,relpred);
dc = c(2)-c(1);
Xtot = sum(X)*dc;
fit_fixed_time = (1-alpha1*c)*(1-Xtot) - (1-alpha2*c)*Y./(ph+xbar) - m1;
f3 = figure;
plot(c,fit_fixed_time,'--o')
xlabel("Prey trait")
ylabel("Prey fitness")
% save("num_sim_data/growth_rate_vs_c_data_alpha1"+strrep(num2str(alpha1),".","dot")+"_alpha2"+strrep(num2str(alpha2),".","dot")+"_m1"+strrep(num2str(m1),".","dot")...
%     +"_m2"+strrep(num2str(m2),".","dot"),"c","fit","alpha1","alpha2","m1","m2","d","ph","gamma", "v_op", "t")


%% plot all together
% try
% alpha1_col = [0.2,0.29,0.38];
% f4 = figure;
% 
% for aa = 1:length(alpha1_col)
%     alpha1 = alpha1_col(aa);
%     load("num_sim_data/growth_rate_vs_c_data_alpha1"+strrep(num2str(alpha1),".","dot")+"_alpha2"+strrep(num2str(alpha2),".","dot")+"_m1"+strrep(num2str(m1),".","dot")...
%         +"_m2"+strrep(num2str(m2),".","dot"))
%     subplot(1,2,1)
%     hold on 
%     grid on
%     p(aa) = plot(c,fit,'-', "DisplayName","$\alpha_1 = "+num2str(alpha1)+"$");
%     subplot(1,2,2)
%     hold on
%     grid on
%     plotind_single = t(end); % end point of simulation
%     plot(c,v_op(1:M)/((c(2)-c(1))*sum(v_op(1:M))))
% end
% subplot(1,2,1)
% leg = legend(p,"Location","Southeast");
% set(leg,"Interpreter","latex")
% xlabel("Trait value, $c$", "Interpreter","latex")
% ylabel("Individual prey growth rate, $g(X,Y;c)/X$", "Interpreter","latex")
% subplot(1,2,2)
% xlabel("Trait value, c")
% ylabel("Prey frequency")
% %sgtitle("$m_1 = " + num2str(m1) +", m_2 = "+num2str(m2)+", \alpha_2 = "+num2str(alpha2)+  ", p = " + num2str(ph) + ", d = "+ num2str(d)   + "$", 'interpreter', 'latex')
% set(f4,'Windowstyle','normal')
% set(findall(f4,'-property','FontSize'),'FontSize',11)
% set(f4,'Units','centimeters')
% set(f4,'Position',[18 1 18 9])
% end

%% fitness figures
if const == 0
    f5 = figure(15);
    X = v_op(:,1:M);
    Y = v_op(:,end);
    relpred = (1-alpha2*c).*X;
    xbar = NaN*ones(1,height(X)); Xtot = xbar;
    dc = c(2)-c(1);
    for xx = 1:height(X)
        xbar(xx) = trapz(c,relpred(xx,:));
        Xtot(xx) = sum(X(xx,:))*dc;
        fit(xx,:) = (1-alpha1*c)*(1-Xtot(xx)) - (1-alpha2*c)*Y(xx)./(ph+xbar(xx)) - m1;
    end
%     maxfit = max(fit(:)); minfit = min(fit(:));
%     for xx = 1:height(X)
%         fit(xx,:) = (fit(xx,:)-minfit)/(maxfit-minfit);
%     end
    eqfit_ind = find(islocalmin(std(fit')));
    
    contourf(c,t_op,fit, 'linestyle', 'none')
    colormap(flipud(summer))
    shading interp
    colorbar
    hold on
    for aa = 1:length(eqfit_ind)
        yline(t_op(eqfit_ind(aa)), '--r','LineWidth',2)
    end
    ylabel("Time")
    xlabel("Prey defence trait")
    title("$\alpha_1 = "+num2str(alpha1)+"$","interpreter","latex")


    set(f5,'Windowstyle','normal')
    set(findall(f5,'-property','FontSize'),'FontSize',11)
    set(f5,'Units','centimeters')
    set(f5,'Position',[18 1 17/2 17/2])
    figname = "../../Ecol_paper/figures/DD"+strrep("_d"+num2str(d)+"_ph"+num2str(ph)+"_gamma"+num2str(gamma)+...
        "_alpha1"+num2str(alpha1)+"_alpha2"+num2str(alpha2)+"_m1"+num2str(m1)+"_m2"+num2str(m2),'.','dot')+".eps";

    % exportgraphics(f5,figname,"Resolution",2000)

end
