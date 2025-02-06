%% Prey defence Hopf locus test 
% This script uses direct numerical simulations to calculate the location of
% the transition between stable and oscillatory solutions.

clear; 
close all;
plotonly = 1;
para1 = 'alpha2';
para2 = 'alpha1';
alt = 0;
if alt == 1
    alttext = "_alt";
elseif alt == 2
    alttext = "_alt_superlinear_cost";
elseif alt == 3
    alttext = "_alt_sat_eff";
else
    alttext = "";
end
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
s = 0.5; % saturation level of prey defence efficiency

%% Mesh

tmax = 1000; %Integration range for solver
M = 2^8; %Number of trait points





filename = "num_sim_data/hopf_sim_data_"+para1+ "_" + para2 + "_change"+strrep("_d"+num2str(d)+"_ph"+num2str(ph)+"_gamma"+num2str(gamma)+...
        "_alpha1"+num2str(alpha1)+"_alpha2"+num2str(alpha2)+"_m1"+num2str(m1)+"_m2"+num2str(m2)+alttext,'.','dot');

if plotonly ~=1
%% ODE Solver
    options = odeset('Stats', 'off','MaxStep',1e-2,'NonNegative',1:M+1); %Note v = [u;w] = [plants; water]
  
    para2_col_alt = [0.1,0.2,0.25];
%     para2_col_alt(mod(100*para2_col_alt,5)==0) = [];
    for bb = 1:length(para2_col_alt)
    disp("Major step "+num2str(bb)+" of "+num2str(length(para2_col_alt)))
    para2_col = para2_col_alt(bb);
    para1_col = linspace( 5*para2_col-0.1,5*para2_col+0.1,21);
%     para1_col = 0.1:0.1:2.0;
%     para2_col = 0.1:0.1:2;

    limitcycle = zeros(length(para1_col),length(para2_col));
    for aa = 1:length(para1_col)
        eval([para1 '=para1_col(aa);'])
        for mm = 1:length(para2_col)
            eval([para2 '=para2_col(mm);'])
            disp("Step "+num2str((aa-1)*length(para2_col)+mm)+ " of "+ num2str(length(para1_col)*length(para2_col)));
            if alt == 2
                cmax = sqrt(1/alpha1);
            else
                cmax = 1/alpha1; %Space domain size of interest
            end
            c=linspace(0,cmax,M);         
            %% IC
            u0 = 0.5*ones(1,length(c));%0.1+0.1*rand(1,length(c));
            u0(end+1) = 0.1;
            try
            [t,v,totalprey,medianc,meanc,L,v_op,totalprey_op,t_op,medianc_op,meanc_op,interq_trait_op,phaselag_prey_pred,phaselag_pred_trait] = prey_defence_single_run_fun(c,M,d,alpha1,alpha2,ph,gamma,m2,m1,tmax,u0,options,alt,s);
    
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
        save(filename, 'datacol')    
    end
% save(filename, 'datacol')    
    
else
    load(filename);
end
%%
fields = fieldnames(datacol);
% rmind = [];
for ff = 1:length(fields)
   data = datacol.(fields{ff});
%    if any(abs(data.para1 - 0.35)<1e-3)
%        rmind = [rmind,ff];
%    end
   [para1_col,para2_col] = meshgrid(datacol.(fields{ff}).para1, datacol.(fields{ff}).para2);
   if size(para2_col) ~= size(para1_col)
       para2_col = repmat(para2_col, width(para1_col),1);
       para2_col = para2_col';
   end
   limitcycle = datacol.(fields{ff}).limitcycle';
   plot(para2_col(limitcycle==1), para1_col(limitcycle==1), 'o', 'color', 'r')
   plot(para2_col(limitcycle==0), para1_col(limitcycle==0), 'o', 'color', 'k')
end

grid on
ylabel(para1)
xlabel(para2)
save(filename, 'datacol')

%% plot boundary only
f1 = figure;
ms = 3;
fields = fieldnames(datacol);
para1_b = []; para2_b = [];
para1_all = []; para2_all = []; limitcycle_all = [];
for ff = 1:length(fields)
   data = datacol.(fields{ff});
   [P1,P2] = meshgrid(data.para1,data.para2);
   para1_all = [para1_all, P1(:)']; para2_all = [para2_all, P2(:)'];
   lc = data.limitcycle';
   limitcycle_all = [limitcycle_all, lc(:)'];

%    para1_col = repmat(datacol.(fields{ff}).para1, width(datacol.(fields{ff}).para2),1);
%    para2_col  = datacol.(fields{ff}).para2;
%    if size(para2_col) ~= size(para1_col)
%        para2_col = repmat(para2_col, width(para1_col),1);
%        para2_col = para2_col';
%    end
%    limitcycle = datacol.(fields{ff}).limitcycle';
%    for aa = 1:length(data.para1) 
%         b_ind = find(diff(data.limitcycle(aa,:))~=0);
%         
%         for bb = 1:length(b_ind)
%             para1_b = [para1_b, data.para1(aa)];
%             para2_b = [para2_b, mean(data.para2(b_ind(bb):b_ind(bb)+1))];
%         end      
%         
%    end
end

para1_plot = linspace(0.1,2,1e3); para2_plot = linspace(0.1,2,1e3);
[P1_plot,P2_plot] = meshgrid(para1_plot,para2_plot);
limitcycle_plot = NaN*ones(1,length(P1_plot(:)));
for pp = 1:length(P1_plot(:))
    [~,minind] = min(abs((P1_plot(pp)-para1_all).^2 + (P2_plot(pp)-para2_all).^2));
    limitcycle_plot(pp) = limitcycle_all(minind);
end

% [para1_b,sortind] = sort(para1_b); para2_b = para2_b(sortind);
% [para1_b,ia,ic] = unique(para1_b);

% plot(P2_plot(limitcycle_plot==1), P1_plot(limitcycle_plot==1), 'o', 'color', 'r')
% hold on
% plot(P2_plot(limitcycle_plot==0), P1_plot(limitcycle_plot==0), 'o', 'color', 'k')

contour(P2_plot,P1_plot,reshape(limitcycle_plot,size(P1_plot)),1,"-k", "LineWidth",2)

% plot(para2_b(ia),para1_b,'o', 'color', 'k', 'MarkerSize',ms)

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
% if strcmp(para1,'m2')
%     xlim([0,1.5])
% end
if alt == 0
    xlim([0.1,1])
    ylim([0.1,1])
else
    xlim([0.1,2])
    ylim([0.1,2])
end
pbaspect([1 1 1])

set(f1,'Windowstyle','normal')
% set(f1,'Resize','off')
set(findall(f1,'-property','FontSize'),'FontSize',11)
set(f1,'Units','centimeters')
set(f1,'Position',[18 1 5.5 5.5])

% saveas(f1,"Ecol_paper/figures/hopf_locus_num_sim_"+para1+"_"+para2, 'epsc')

hold on
xplot = linspace(0,2);
% plot(xplot,0.8*xplot,'--','LineWidth',1)
% plot(xplot,2*xplot,'--','LineWidth',1)
% plot(xplot,0.15+0.7*xplot,'--','LineWidth',1)


saveas(f1,"../../Ecol_paper/figures/hopf_locus_num_sim_"+para1+"_"+para2+"_trade_offs"+alttext, 'epsc')