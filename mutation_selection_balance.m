%% Prey defence mutation-selection-balance calculation
% This script calulates the individual variability due to a
% mutation-selection balance for a wide range of parameters. 

clear; close all;

alt = 1;
if alt == 1
    alttext = "_alt";
else
    alttext = "";
end

%% Parameters
d=0.001; dtemp = d;% mutation rate
alpha1 = 0.75; alpha1temp = alpha1; % cost of prey defence
alpha2 = 0.5; alpha2temp = alpha2; % prey defence efficiency
m1 = 0.2; %prey mortality
m2 = 0.2; m2temp = m2; %pred mortality (LV only)
ph = 0.5; %predation half saturation constant (extension only)
gamma = 4; % prey to predator conversion
s =0.5;
switchparachoice = "d";
% switchpara = [logspace(-5,-4,10),logspace(-4,-3,10),logspace(-3,-2,10),logspace(-2,-1,10)];
if switchparachoice == "m2"
    switchpara = linspace(0,3,30);
elseif switchparachoice == "p"
    switchpara = 0:0.01:1; 
elseif switchparachoice == "d"
    switchpara = logspace(-6,-1,100);
elseif switchparachoice == "gamma"
    switchpara = 0:0.05:8; 
elseif switchparachoice == "alpha1"
    switchpara = 0.05:0.05:2;
else
    switchpara = 0.05:0.05:2;
end
%% Mesh
tmax = 1000; %Integration range for solver
M = 2^8; %Number of trait points





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
    cmax = min([1/alpha1]); %Space domain size of interest
    c=linspace(0,cmax,M);
    % IC
    u0 = 0.5*ones(1,length(c));%0.1+0.1*rand(1,length(c));
    u0(end+1) = 0.0; % no predators in this case
    [t,v,totalprey,medianc,meanc,L,v_op,totalprey_op,t_op,medianc_op,meanc_op,varc_trait_op,phaselag_prey_pred,phaselag_pred_trait,phaselag_mean_var,varc] = prey_defence_single_run_fun(c,M,d,alpha1,alpha2,ph,gamma,m2,m1,tmax,u0,options,alt,s);



    %% calc IQR
%     np = length(c);
%         interq_trait=zeros(1,length(t)); traitdist = zeros(np,length(t));
%         for i=1:length(t) % loop through all times
%             traitdist(1:np,i) = v(i,1:np)/sum(v(i,1:np)); % frequency of trait value
%             csum = cumsum(traitdist(:,i)); %cdf of trait dist
%             cspace = linspace(0,1,1000); % larger c vector
%             csum = interp1(c',csum,cspace); % interpolate onto larger c vector
%             p25ind = find(csum<0.25); % find 25 percentile
%             if ~isempty(p25ind) 
%                 p25ind = p25ind(end)+1; cp25 = cspace(p25ind);
%             else
%                 cp25 = cspace(1);
%             end
%             p75ind = find(csum<0.75); % find 75 percentile
%             if ~isempty(p75ind)
%                 p75ind = p75ind(end); cp75 = cspace(p75ind);
%             else
%                 cp75 = cspace(end);
%             end
%             interq_trait(i) = cp75-cp25; % interquartile range of traits
%         end
        mut_sel_varc(ss) = varc(end);
end
switch switchparachoice
    case "m2"
        m2 = m2temp;
    case "alpha1"
        alpha1 = alpha1temp;
    case "alpha2"
        alpha2 = alpha2temp;
    case "d"
        d = dtemp;
end


filename = "num_sim_data/mut_sel_balance_data_"+switchparachoice+"_change"+strrep("_d"+num2str(d)+"_ph"+num2str(ph)+"_gamma"+num2str(gamma)+...
        "_alpha1"+num2str(alpha1)+"_alpha2"+num2str(alpha2)+"_m1"+num2str(m1)+"_m2"+num2str(m2)+alttext,'.','dot');

save(filename,"switchpara","mut_sel_varc")

if switchparachoice == "d"
f = figure;
semilogx(switchpara,mut_sel_varc, '--o')
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