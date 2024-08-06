%% Prey defence eco-evol model: Bifurcation diagrams
% This script calculates and visualises the bifurcation diagrams.

clear; 
close all;
plotonly = 1;

%% set which parameter is changed

parachange = "alpha2";
if parachange == "m2"
    para_col = 0.06:0.02:3;
elseif parachange == "d"
    para_col = logspace(-6,-2,40);
else
    para_col = 0:0.01:1; 
end
%% Parameters
d=0.001; % mutation rate
alpha1 = 0.75; % max growth
alpha2 = 0.5; % max predation (LV) or max predation = pmax*alpha2 (extension)
m1 = 0.2; %prey mortality
m2 = 0.2; %pred mortality (LV only)
ph = 0.5; %predation half saturation constant (extension only)
gamma = 4; % prey to predator conversion
filename = "num_sim_data/"+parachange+"_change"+strrep("_d"+num2str(d)+"_ph"+num2str(ph)+"_gamma"+num2str(gamma)+...
        "_alpha1"+num2str(alpha1)+"_alpha2"+num2str(alpha2)+"_m1"+num2str(m1)+"_m2"+num2str(m2),'.','dot');
filename1 = "num_sim_data/mut_sel_balance_data_"+parachange+"_change"+strrep("_d"+num2str(d)+"_ph"+num2str(ph)+"_gamma"+num2str(gamma)+...
        "_alpha1"+num2str(alpha1)+"_alpha2"+num2str(alpha2)+"_m1"+num2str(m1)+"_m2"+num2str(m2),'.','dot');

if plotonly == 0
    %% Mesh
    cmax = 1; %Space domain size of interest
    tmax = 1000; %Integration range for solver
    M = 2^8; %Number of trait points
    if M>1
        c=linspace(0,cmax,M);
    else
        c=0.0;
    end
    
    %% IC
    % u0 = [100*rand(1,length(c))/length(c),10*rand];
    u0 = 0.5*ones(1,length(c));%0.1+0.1*rand(1,length(c));
    % if M>1 
    %     u0(c<1/3) = 0; u0(c>2/3) = 0;
    % end
    u0(end+1) = 0.5;
    
    options = odeset('Stats', 'off','MaxStep',1e-2,'NonNegative',1:M+1); 
    
    
    meanprey = NaN*ones(1,length(para_col)); minprey = meanprey; maxprey = meanprey; meanpred = meanprey; minpred = meanprey; maxpred = meanprey; meanmeantrait = meanprey;
    minmeantrait = meanprey; maxmeantrait = meanprey; meanvartrait = meanprey; minvartrait = meanprey; maxvartrait = meanprey; wavelength = meanprey;
    lag_prey_pred = meanprey; lag_pred_trait = meanprey;
    for mm = 1:length(para_col)
        fprintf("\n Step " + num2str(mm) + " of "+ num2str(length(para_col)))
        if parachange == "m2"
            m2 = para_col(mm);
        elseif parachange == "alpha2"
            alpha2 = para_col(mm);
        elseif parachange == "alpha1"
            alpha1 = para_col(mm);
        elseif parachange == "m1"
            m1 = para_col(mm);
        elseif parachange == "d"
            d = para_col(mm);
        else
            error("Not a valid parameter to change")
        end
        [t,v,totalprey,medianc,meanc,L,v_op,totalprey_op,t_op,medianc_op,meanc_op,varc_op,phaselag_prey_pred,phaselag_pred_trait,phaselag_mean_iqr] = prey_defence_single_run_fun(c,M,d,alpha1,alpha2,ph,gamma,m2,m1,tmax,u0,options);
        
        meanprey(mm) = mean(totalprey_op); minprey(mm) = min(totalprey_op); maxprey(mm) = max(totalprey_op);
        meanpred(mm) = mean(v_op(:,end)); minpred(mm) = min(v_op(:,end)); maxpred(mm) = max(v_op(:,end));
        meanmeantrait(mm) = mean(meanc_op); minmeantrait(mm) = min(meanc_op); maxmeantrait(mm) = max(meanc_op);
        meanvartrait(mm) = mean(varc_op); minvartrait(mm) = min(varc_op); maxvartrait(mm) = max(varc_op);
        wavelength(mm) = L;
        lag_prey_pred(mm) = phaselag_prey_pred; lag_pred_trait(mm) = phaselag_pred_trait; lag_mean_iqr(mm) = phaselag_mean_iqr;
    
    end
    
    save(filename,"para_col","meanprey","maxprey","minprey","meanpred","minpred","maxpred","meanmeantrait","minmeantrait","maxmeantrait",...
        "meanvartrait","minvartrait","maxvartrait","wavelength","lag_pred_trait","lag_prey_pred","lag_mean_iqr")
else
    load(filename)
end
col = lines;
meancol = 'k';
maxcol = col(1,:);
mincol = col(2,:);
ms = 3; % markersize
f = figure;

load(filename1)



subplot(2,2,1)
hold on
grid on
plot(para_col,minprey,'--o','color', mincol, 'MarkerSize',ms)
plot(para_col,maxprey,'--o','color', maxcol, 'MarkerSize',ms)
plot(para_col,meanprey,'--o','color', meancol, 'MarkerSize',ms)
ylabel("Prey biomass")
xlim([para_col(1),para_col(end)])
ylim([0,1])
pbaspect([1.5 1 1])
% title("A")
% ax=gca;
% ax.TitleHorizontalAlignment = 'left'; 
if parachange == "d"
    set(gca,'xscale','log')
end

subplot(2,2,2)
hold on
grid on
plot(para_col,minpred,'--o','color', mincol, 'MarkerSize',ms)
plot(para_col,maxpred,'--o','color', maxcol, 'MarkerSize',ms)
plot(para_col,meanpred,'--o','color', meancol, 'MarkerSize',ms)
ylabel("Pred biomass")
xlim([para_col(1),para_col(end)])
ylim([0,2])
pbaspect([1.5 1 1])
% title("B")
% ax=gca;
% ax.TitleHorizontalAlignment = 'left'; 
if parachange == "d"
    set(gca,'xscale','log')
end

subplot(2,2,3)
hold on
grid on
plot(para_col,minmeantrait,'--o','color', mincol, 'MarkerSize',ms)
plot(para_col,maxmeantrait,'--o','color', maxcol, 'MarkerSize',ms)
plot(para_col,meanmeantrait,'--o','color', meancol, 'MarkerSize',ms)
ylabel("Mean trait")
xlim([para_col(1),para_col(end)])
ylim([0,1])
pbaspect([1.5 1 1])
% title("C")
% ax=gca;
% ax.TitleHorizontalAlignment = 'left'; 
if parachange == "d"
    set(gca,'xscale','log')
end

subplot(2,2,4)
hold on
grid on
plot(para_col,sqrt(minvartrait),'--o','color', mincol, 'MarkerSize',ms)
plot(para_col,sqrt(maxvartrait),'--o','color', maxcol, 'MarkerSize',ms)
plot(para_col,sqrt(meanvartrait),'--o','color', meancol, 'MarkerSize',ms)
plot(switchpara,sqrt(mut_sel_varc),'--', 'color', col(5,:))
ylabel("Std")
xlim([para_col(1),para_col(end)])
ylim([0,0.55])
pbaspect([1.5 1 1])
% title("D")
% ax=gca;
% ax.TitleHorizontalAlignment = 'left';
if parachange == "d"
    set(gca,'xscale','log')
end
if parachange == "m2"
    xlabel("Predator mortality, $m_2$", "Interpreter","latex", "Position", [-0.25  -0.1])
elseif parachange == "alpha2"
    xlabel("Prey defence efficiency, $\alpha_2$", "Interpreter","latex", "Position", [-0.15  -0.1])
elseif parachange == "alpha1"
    xlabel("Prey defence cost, $\alpha_1$", "Interpreter","latex", "Position", [-0.15  -0.1])
elseif parachange == "m1"
    xlabel("Prey mortality, $m_1$", "Interpreter","latex", "Position", [-0.15 -0.1])
elseif parachange == "d"
    xlabel("Prey mutation strength, $d$", "Interpreter","latex", "Position", [-0.15 -0.1])
    sgtitle("$\alpha_1 = "+num2str(alpha1)+"$","Interpreter","latex")
end

set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[18 1 13 9])

% saveas(f,"../../Ecol_paper/figures/bif_diag_"+parachange+"_numsim", 'epsc')


f1 = figure;
subplot(1,2,1)
hold on
grid on
plot(para_col,wavelength,'--o','color', meancol, 'MarkerSize',ms)
ylabel("Wavelength")
xlim([para_col(1),para_col(end)])
ylim([0,100])
% title("E")
% ax=gca;
% ax.TitleHorizontalAlignment = 'left'; 
pbaspect([1.5 1 1])
if parachange == "d"
    set(gca,'xscale','log')
end


subplot(1,2,2)
hold on
grid on

lag_prey_pred(lag_prey_pred > 0.5) = 1-lag_prey_pred(lag_prey_pred>0.5);
plot(para_col,lag_prey_pred,'--o','color', meancol, 'MarkerSize',ms)
ylabel("Phase lag")
xlim([para_col(1),para_col(end)])
ylim([0,0.5])
pbaspect([1.5 1 1])

lag_pred_trait(lag_pred_trait > 0.5) = 1-lag_pred_trait(lag_pred_trait>0.5);
plot(para_col,lag_pred_trait,'--o','color', col(1,:), 'MarkerSize',ms)

lag_mean_iqr(lag_mean_iqr > 0.5) = 1-lag_mean_iqr(lag_mean_iqr>0.5);
plot(para_col,lag_mean_iqr,'--o','color', col(2,:), 'MarkerSize',ms)

if parachange == "d"
    set(gca,'xscale','log')
end

xlabel("Test")
if parachange == "m2"
    xlabel("Predator mortality, $m_2$", "Interpreter","latex", "Position", [-0.25  -0.08])
elseif parachange == "alpha2"
    xlabel("Prey defence efficiency, $\alpha_2$", "Interpreter","latex", "Position", [-0.15  -0.08])
elseif parachange == "alpha1"
    xlabel("Prey defence cost, $\alpha_1$", "Interpreter","latex", "Position", [-0.15  -0.08])
elseif parachange == "m1"
    xlabel("Prey mortality, $m_1$", "Interpreter","latex", "Position", [-0.15 -0.08])
elseif parachange == "d"
    xlabel("Prey mutation strength, $d$", "Interpreter","latex", "Position", [-0.15 -0.1])
end


set(f1,'Windowstyle','normal')
set(findall(f1,'-property','FontSize'),'FontSize',11)
set(f1,'Units','centimeters')
set(f1,'Position',[18 1 13 4.5])

% saveas(f1,"Ecol_paper/figures/bif_diag_"+parachange+"_numsim_suppl", 'epsc')