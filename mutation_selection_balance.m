%% Prey defence mutation-selection-balance calculation
% This script calulates the individual variability due to a
% mutation-selection balance for a wide range of parameters. 

clear; close all;

%% Parameters
d=0.001; % mutation rate
alpha1 = 0.5; % cost of prey defence
alpha2 = 0.5; % prey defence efficiency
m1 = 0.2; %prey mortality
m2 = 0.5; %pred mortality (LV only)
ph = 0.5; %predation half saturation constant (extension only)
gamma = 4; % prey to predator conversion
switchparachoice = "d";
% switchpara = [logspace(-5,-4,10),logspace(-4,-3,10),logspace(-3,-2,10),logspace(-2,-1,10)];
switchpara = linspace(0,1,10);
%% Mesh
cmax = 1; %Space domain size of interest
tmax = 1000; %Integration range for solver
M = 2^8; %Number of trait points
c=linspace(0,cmax,M);


%% IC
u0 = 0.5*ones(1,length(c));%0.1+0.1*rand(1,length(c));
u0(end+1) = 0.0; % no predators in this case


%% ODE Solver
options = odeset('Stats', 'off','MaxStep',1e-2,'NonNegative',1:M+1); 

for ss = 1:length(switchpara)
    disp("Step "+num2str(ss)+" of "+num2str(length(switchpara)))
    switch switchparachoice
        case "m2"
            m2 = switchpara(ss);
        case "alpha1"
            alpha1 = switchpara(ss);
        case "alpha2"
            alpha2 = switchpara(ss);
        case "d"
            d = switchpara(ss);
    end
    [t,v,totalprey,medianc,meanc,L,v_op,totalprey_op,t_op,medianc_op,meanc_op,interq_trait_op,phaselag_prey_pred,phaselag_pred_trait] = prey_defence_single_run_fun(c,M,d,alpha1,alpha2,ph,gamma,m2,m1,tmax,u0,options);



    %% calc IQR
    np = length(c);
        interq_trait=zeros(1,length(t)); traitdist = zeros(np,length(t));
        for i=1:length(t) % loop through all times
            traitdist(1:np,i) = v(i,1:np)/sum(v(i,1:np)); % frequency of trait value
            csum = cumsum(traitdist(:,i)); %cdf of trait dist
            cspace = linspace(0,1,1000); % larger c vector
            csum = interp1(c',csum,cspace); % interpolate onto larger c vector
            p25ind = find(csum<0.25); % find 25 percentile
            if ~isempty(p25ind) 
                p25ind = p25ind(end)+1; cp25 = cspace(p25ind);
            else
                cp25 = cspace(1);
            end
            p75ind = find(csum<0.75); % find 75 percentile
            if ~isempty(p75ind)
                p75ind = p75ind(end); cp75 = cspace(p75ind);
            else
                cp75 = cspace(end);
            end
            interq_trait(i) = cp75-cp25; % interquartile range of traits
        end
        mut_sel_iqr(ss) = interq_trait(end);
end

switch switchparachoice
    case "m2"
        filename = "num_sim_data/mut_sel_balance_data_m2";
    case "alpha1"
        filename = "num_sim_data/mut_sel_balance_data_alpha1";
    case "alpha2"
        filename = "num_sim_data/mut_sel_balance_data_alpha2";
    case "d"
        filename = "num_sim_data/mut_sel_balance_data_d";
end
save(filename,"switchpara","mut_sel_iqr")

if switchparachoice == "d"
f = figure;
semilogx(switchpara,mut_sel_iqr, '--o')
hold on
grid on
xlabel("Trait difusion, $d$", "Interpreter", "latex")
ylabel("IQR due to mutation-selection balance", "Interpreter", "latex")
set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[18 1 18 9])
xlim([1e-5,1e-1])
end