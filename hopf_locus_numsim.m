%% Prey defence Hopf locus test 
% This script uses direct numerical simulations to calculate the location of
% the transition between stable and oscillatory solutions.

clear; 
close all;
plotonly = 1;
para1 = 'm2';
para2 = 'alpha2';
f = figure;
hold on
%% Parameters
d=0.001; % mutation rate
alpha1 = 0.75; % max growth
m1 = 0.2; %prey mortality
m2 = 0.2;
alpha2 = 0.5;
ph = 0.5; %predation half saturation constant (extension only)
gamma = 4; % prey to predator conversion

%% Mesh
cmax = 1; %Space domain size of interest
tmax = 1000; %Integration range for solver
M = 2^8; %Number of trait points
c=linspace(0,cmax,M);


%% IC
u0 = 0.5*ones(1,length(c));%0.1+0.1*rand(1,length(c));
u0(end+1) = 0.5;

filename = "num_sim_data/hopf_sim_data_"+para1+ "_" + para2 + "_change"+strrep("_d"+num2str(d)+"_ph"+num2str(ph)+"_gamma"+num2str(gamma)+...
        "_alpha1"+num2str(alpha1)+"_alpha2"+num2str(alpha2)+"_m1"+num2str(m1)+"_m2"+num2str(m2),'.','dot');

if plotonly ~=1
%% ODE Solver
    options = odeset('Stats', 'off','MaxStep',1e-2,'NonNegative',1:M+1); %Note v = [u;w] = [plants; water]
  
    para1_col = 0.4:0.05:0.9;
    para2_col = 0.025:0.05:0.875;

    limitcycle = zeros(length(para1_col),length(para2_col));
    for aa = 1:length(para1_col)
        eval([para1 '=para1_col(aa);'])
        for mm = 1:length(para2_col)
            eval([para2 '=para2_col(mm);'])
            disp("Step "+num2str((aa-1)*length(para2_col)+mm)+ " of "+ num2str(length(para1_col)*length(para2_col)));
            try
            [t,v,totalprey,medianc,meanc,L,v_op,totalprey_op,t_op,medianc_op,meanc_op,interq_trait_op,phaselag_prey_pred,phaselag_pred_trait] = prey_defence_single_run_fun(c,M,d,alpha1,alpha2,ph,gamma,m2,m1,tmax,u0,options);
    
            if max(abs(max(totalprey_op)-min(totalprey_op))./mean(totalprey_op)) > 1e-2 && max(totalprey_op) > 1e-2
                disp("Limit cycle solution detected /n")
                limitcycle(aa,mm) = 1;
            end
            end
        end
    end
        
    
    try
        load(filename);
        fields = fieldnames(datacol);
        datacol.("Data_"+num2str(length(fields)+1)).para1 = para1_col;
        datacol.("Data_"+num2str(length(fields)+1)).para2 = para2_col;
        datacol.("Data_"+num2str(length(fields)+1)).limitcycle = limitcycle;
    catch
        disp("First run")
        datacol.("Data_1").para1 = para1_col;
        datacol.("Data_1").para2 = para2_col;
        datacol.("Data_1").limitcycle = limitcycle;
    end

    
    
else
    load(filename);
end
%%
fields = fieldnames(datacol);
for ff = 1:length(fields)
   data = datacol.(fields{ff});
   [para1_col,para2_col] = meshgrid(datacol.(fields{ff}).para1, datacol.(fields{ff}).para2);
   if size(para2_col) ~= size(para1_col)
       para2_col = repmat(para2_col, width(para1_col),1);
       para2_col = para2_col';
   end
   limitcycle = datacol.(fields{ff}).limitcycle';
   plot(para1_col(limitcycle==1),para2_col(limitcycle==1), 'o', 'color', 'r')
   plot(para1_col(limitcycle==0),para2_col(limitcycle==0), 'o', 'color', 'k')
end

grid on
xlabel(para1)
ylabel(para2)
save(filename, 'datacol')

%% plot boundary only
f1 = figure;
ms = 3;
fields = fieldnames(datacol);
para1_b = []; para2_b = [];
for ff = 1:length(fields)
   data = datacol.(fields{ff});
   para1_col = repmat(datacol.(fields{ff}).para1, width(datacol.(fields{ff}).para2),1);
   para2_col  = datacol.(fields{ff}).para2;
   if size(para2_col) ~= size(para1_col)
       para2_col = repmat(para2_col, width(para1_col),1);
       para2_col = para2_col';
   end
   limitcycle = datacol.(fields{ff}).limitcycle';
   for aa = 1:length(data.para1) 
        b_ind = find(diff(data.limitcycle(aa,:))~=0);
        
        for bb = 1:length(b_ind)
            para1_b = [para1_b, data.para1(aa)];
            para2_b = [para2_b, mean(data.para2(b_ind(bb):b_ind(bb)+1))];
        end      
        
   end
end
[para1_b,sortind] = sort(para1_b); para2_b = para2_b(sortind);
[para1_b,ia,ic] = unique(para1_b);
plot(para2_b,para1_b,'o', 'color', 'k', 'MarkerSize',ms)

if strcmp(para1,'m2') && strcmp(para2,'m1')
    hold on
    m1plot = linspace(0,3);
    m2plot = gamma*(1-ph-m1plot)./(1+ph+m1plot);
    plot(m2plot,m1plot)
    legend("Eco-evol.", "Eco only")
end
grid on
if strcmp(para1,'alpha2')
    ylabel("Prey defence efficiency, $\alpha_2$", "Interpreter","latex")
elseif strcmp(para1,'m2')
    ylabel("Predator mortality, $m_2$", "Interpreter","latex")
else
    ylabel(para1)
end

if strcmp(para2,'alpha1')
    xlabel("Prey defence cost, $\alpha_1$", "Interpreter","latex")
elseif strcmp(para2,'m1')
    xlabel("Prey mortality, $m_1$", "Interpreter","latex")
elseif strcmp(para2,'alpha2')
    xlabel("Prey defence efficiency, $\alpha_2$", "Interpreter","latex")
else
    xlabel(para2)
end
xlim([f.Children.XLim])
ylim([f.Children.YLim])
if strcmp(para1,'m2')
    xlim([0,1.5])
end
pbaspect([1 1 1])

set(f1,'Windowstyle','normal')
% set(f1,'Resize','off')
set(findall(f1,'-property','FontSize'),'FontSize',11)
set(f1,'Units','centimeters')
set(f1,'Position',[18 1 8.5 8.5])

% saveas(f1,"Ecol_paper/figures/hopf_locus_num_sim_"+para1+"_"+para2, 'epsc')

hold on
xplot = linspace(0,1);
plot(xplot,xplot,'--')
plot(xplot,2*xplot,'--')
plot(xplot,0.25+0.7*xplot,'--')

xlim([0,1])
ylim([0,1])
% saveas(f1,"Ecol_paper/figures/hopf_locus_num_sim_"+para1+"_"+para2+"_trade_offs", 'epsc')