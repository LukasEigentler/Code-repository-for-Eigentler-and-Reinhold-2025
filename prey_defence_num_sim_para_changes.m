%% Prey defence eco-evol model: Bifurcation diagrams
% This script calculates and visualises the bifurcation diagrams.

clear; 
close all;
plotonly = 1;
alt = 1; % 0 for model in which f1 multiplies the net growth term; 1 for model in which f1 multiplies growth rate only; 2 for superlinear cost; 3 for saturation of efficiency at s
if alt == 1
    alttext = "_alt";
else
    alttext = "";
end
%% set which parameter is changed

parachange = "p";
if parachange == "m2"
    para_col = 0.1:0.02:3;

elseif parachange == "alpha1"
    if alt == 1
%         para_col =[0.1:0.01:0.47, 0.472:0.002:0.49, 0.5:0.01:2];
        para_col = [0.1:0.002:0.4, 0.41:0.01:2];
    else
        para_col =0.1:0.01:2;
    end
elseif parachange == "alpha2"
    if alt == 1
        para_col = [0.01:0.01:0.77, 0.772:0.002:0.78, 0.79:0.01:2];
    else
        para_col =0.01:0.01:2;
    end
elseif parachange == "p"
    para_col = 0.01:0.01:1; 
elseif parachange == "d"
    para_col = logspace(-6,-1,100);
elseif parachange == "gamma"
    para_col = 0.1:0.05:5; 
end
%% Parameters
d=0.001; % mutation rate
alpha1 = 0.75; % max growth
alpha2 = 0.5; % max predation (LV) or max predation = pmax*alpha2 (extension)
m1 = 0.2; %prey mortality
m2 = 0.2; %pred mortality (LV only)
ph = 0.5; %predation half saturation constant (extension only)
gamma = 4; % prey to predator conversion
s = 0.5;

filename = "num_sim_data/"+parachange+"_change"+strrep("_d"+num2str(d)+"_ph"+num2str(ph)+"_gamma"+num2str(gamma)+...
        "_alpha1"+num2str(alpha1)+"_alpha2"+num2str(alpha2)+"_m1"+num2str(m1)+"_m2"+num2str(m2)+alttext,'.','dot');
filename1 = "num_sim_data/mut_sel_balance_data_"+parachange+"_change"+strrep("_d"+num2str(d)+"_ph"+num2str(ph)+"_gamma"+num2str(gamma)+...
        "_alpha1"+num2str(alpha1)+"_alpha2"+num2str(alpha2)+"_m1"+num2str(m1)+"_m2"+num2str(m2)+alttext,'.','dot');


if plotonly == 0
    %% Mesh
    tmax = 1000; %Integration range for solver
    M = 2^8; %Number of trait points
    
    
    
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
        elseif parachange == "gamma"
            gamma = para_col(mm);
        elseif parachange == "p"
            ph = para_col(mm);
            
        else
            error("Not a valid parameter to change")
        end
        cmax = min([1/alpha1]); %Space domain size of interest
        if M>1
            c=linspace(0,cmax,M);
        else
            c=0.0;
        end
        % IC

        % u0 = [100*rand(1,length(c))/length(c),10*rand];
        u0 = 0.5*ones(1,length(c))/(c(end)-c(1));
        % if M>1 
        %     u0(c<1/3) = 0; u0(c>2/3) = 0;
        % end
        u0(end+1) = 0.5;

        [t,v,totalprey,medianc,meanc,L,v_op,totalprey_op,t_op,medianc_op,meanc_op,varc_op,phaselag_prey_pred,phaselag_pred_trait,phaselag_mean_iqr] = prey_defence_single_run_fun(c,M,d,alpha1,alpha2,ph,gamma,m2,m1,tmax,u0,options,alt,s);
        
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
lw = 0.5; % Linewidth
f = figure;

load(filename1)

nopreyind = find(maxprey<1e-2);
meanmeantrait(nopreyind) = 0;
minmeantrait(nopreyind) = 0;
maxmeantrait(nopreyind) = 0;
meanvartrait(nopreyind) = 0;
minvartrait(nopreyind) = 0;
maxvartrait(nopreyind) = 0;
wavelength(nopreyind) = 0;
lag_pred_trait(nopreyind) = 0;
lag_prey_pred(nopreyind) = 0;
lag_mean_iqr(nopreyind) = 0;

subplot(2,2,1)
hold on
grid on
plot(para_col,minprey,'o','color', mincol, 'MarkerSize',ms, 'LineWidth',lw)
plot(para_col,maxprey,'o','color', maxcol, 'MarkerSize',ms, 'LineWidth',lw)
plot(para_col,meanprey,'o','color', meancol, 'MarkerSize',ms, 'LineWidth',lw)
ylabel("Prey biomass")
xlim([min(para_col),max(para_col)])
ylim([0,0.8])
pbaspect([1.5 1 1])
% title("A")
% ax=gca;
% ax.TitleHorizontalAlignment = 'left'; 
if parachange == "d"
    set(gca,'xscale','log')
    xticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1])
end

subplot(2,2,2)
hold on
grid on
plot(para_col,minpred,'o','color', mincol, 'MarkerSize',ms, 'LineWidth',lw)
plot(para_col,maxpred,'o','color', maxcol, 'MarkerSize',ms, 'LineWidth',lw)
plot(para_col,meanpred,'o','color', meancol, 'MarkerSize',ms, 'LineWidth',lw)
ylabel("Pred biomass")
xlim([min(para_col),max(para_col)])
ylim([0,3.5])
pbaspect([1.5 1 1])
% title("B")
% ax=gca;
% ax.TitleHorizontalAlignment = 'left'; 
if parachange == "d"
    set(gca,'xscale','log')
    xticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1])
end

subplot(2,2,3)
hold on
grid on
plot(para_col,minmeantrait,'o','color', mincol, 'MarkerSize',ms, 'LineWidth',lw)
plot(para_col,maxmeantrait,'o','color', maxcol, 'MarkerSize',ms, 'LineWidth',lw)
plot(para_col,meanmeantrait,'o','color', meancol, 'MarkerSize',ms, 'LineWidth',lw)
ylabel("Mean trait")
xlim([min(para_col),max(para_col)])
if parachange =="alpha1"
    ylim([0,inf])
else
    ylim([0,1/alpha1])
end
pbaspect([1.5 1 1])
% title("C")
% ax=gca;
% ax.TitleHorizontalAlignment = 'left'; 
if parachange == "d"
    set(gca,'xscale','log')
    xticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1])
end

subplot(2,2,4)
hold on
grid on
plot(para_col,sqrt(minvartrait),'o','color', mincol, 'MarkerSize',ms, 'LineWidth',lw)
plot(para_col,sqrt(maxvartrait),'o','color', maxcol, 'MarkerSize',ms, 'LineWidth',lw)
plot(para_col,sqrt(meanvartrait),'o','color', meancol, 'MarkerSize',ms, 'LineWidth',lw)
plot(switchpara,sqrt(mut_sel_varc),'--', 'color', col(5,:))
ylabel("Std")
xlim([min(para_col),max(para_col)])
if parachange =="alpha1"
    ylim([0,inf])
else
    ylim([0,0.5*1/alpha1])
end
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
    xticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1])
    xlabel("Prey mutation strength, $d$", "Interpreter","latex", "Position", [1.0127e-07 -0.2032 -1.0000])
%     sgtitle("$\alpha_1 = "+num2str(alpha1)+"$","Interpreter","latex")
elseif parachange == "gamma"
    xlabel("Max predator growth rate, $\gamma$", "Interpreter","latex", "Position", [-0.15 -0.1])
elseif parachange == "p"
    xlabel("Predation half saturation constant, $p$", "Interpreter","latex", "Position", [-0.15 -0.1])
end

set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[18 1 13 9])

saveas(f,"../../Ecol_paper/figures/bif_diag_"+parachange+"_numsim"+alttext+"_initial_bounds", 'epsc')
saveas(f,"../../Ecol_paper/figures/bif_diag_"+parachange+"_numsim"+alttext+"_initial_bounds", 'png')


f1 = figure;
subplot(1,2,1)
hold on
grid on
plot(para_col,wavelength,'o','color', meancol, 'MarkerSize',ms, 'LineWidth',lw)
ylabel("Wavelength")
xlim([min(para_col),max(para_col)])
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
plot(para_col,lag_prey_pred,'o','color', 'k', 'MarkerSize',ms, 'LineWidth',lw)
ylabel("Phase lag")
xlim([min(para_col),max(para_col)])
ylim([0,0.5])
pbaspect([1.5 1 1])

lag_pred_trait(lag_pred_trait > 0.5) = 1-lag_pred_trait(lag_pred_trait>0.5);
plot(para_col,lag_pred_trait,'o','color', col(4,:), 'MarkerSize',ms, 'LineWidth',lw)

% lag_mean_iqr(lag_mean_iqr > 0.5) = 1-lag_mean_iqr(lag_mean_iqr>0.5);
% plot(para_col,lag_mean_iqr,'o','color', col(2,:), 'MarkerSize',ms, 'LineWidth',lw)

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
elseif parachange == "gamma"
    xlabel("Max predator growth rate, $\gamma$", "Interpreter","latex", "Position", [-0.15 -0.1])
elseif parachange == "p"
    xlabel("Predation half saturation constant, $p$", "Interpreter","latex", "Position", [-0.15 -0.1])
end


set(f1,'Windowstyle','normal')
set(findall(f1,'-property','FontSize'),'FontSize',11)
set(f1,'Units','centimeters')
set(f1,'Position',[18 1 13 4.5])

saveas(f1,"../../Ecol_paper/figures/bif_diag_"+parachange+"_numsim_suppl"+alttext+"_initial_bounds", 'epsc')
saveas(f1,"../../Ecol_paper/figures/bif_diag_"+parachange+"_numsim_suppl"+alttext+"_initial_bounds", 'png')