%% Change trait bounds for fixed parameters
% This script compares the outcome of simulations for identical parameter
% values but different bounds on the trait, starting from [0,1/alpha1].
clear; close all

plotonly = 1;
plotall = 1; % for plotting all data at once
alt = 1; % 0 for model in which f1 multiplies the net growth term; 1 for model in which f1 multiplies growth rate only; 2 for superlinear cost; 3 for saturation of efficiency at s
if alt == 1
    alttext = "_alt";
else
    alttext = "";
end

if plotall == 0
    %% Parameters
    d=1e-3; % mutation rate
    alpha1 = 0.75; % max growth
    alpha2 = 0.5; % max predation (LV) or max predation = pmax*alpha2 (extension)
    m1 = 0.2; %prey mortality
    m2 = 0.2; %pred mortality (LV only)
    ph = 0.35; %predation half saturation constant (extension only)
    gamma = 4; % prey to predator conversion
    s = 0.5;
    
    bound_deviation = linspace(0,3,31);
    
    
    
    filename = "num_sim_data/trait_bound_test"+strrep("_d"+num2str(d)+"_ph"+num2str(ph)+"_gamma"+num2str(gamma)+...
            "_alpha1"+num2str(alpha1)+"_alpha2"+num2str(alpha2)+"_m1"+num2str(m1)+"_m2"+num2str(m2)+alttext,'.','dot');
    
    if plotonly == 0
    
        %% Mesh
        tmax = 1000; %Integration range for solver
        M = 2^8; %Number of trait points
        
        
        
        options = odeset('Stats', 'off','MaxStep',1e-2,'NonNegative',1:M+1); 
        
        
        meanprey = NaN*ones(1,length(bound_deviation)); minprey = meanprey; maxprey = meanprey; meanpred = meanprey; minpred = meanprey; maxpred = meanprey; meanmeantrait = meanprey;
        minmeantrait = meanprey; maxmeantrait = meanprey; meanvartrait = meanprey; minvartrait = meanprey; maxvartrait = meanprey; wavelength = meanprey;
        lag_prey_pred = meanprey; lag_pred_trait = meanprey;
        for mm = 1:length(bound_deviation)
            fprintf("\n Step " + num2str(mm) + " of "+ num2str(length(bound_deviation)))
            
            z=linspace(0-bound_deviation(mm), 1/alpha1+bound_deviation(mm),M);
            % IC
    
            % u0 = [100*rand(1,length(c))/length(c),10*rand];
            u0 = 0.5*ones(1,length(z));%0.1+0.1*rand(1,length(c));
            % if M>1 
            %     u0(c<1/3) = 0; u0(c>2/3) = 0;
            % end
            u0(end+1) = 0.5;
    
            [t,v,totalprey,medianc,meanc,L,v_op,totalprey_op,t_op,medianc_op,meanc_op,varc_op,phaselag_prey_pred,phaselag_pred_trait,phaselag_mean_iqr] = prey_defence_single_run_fun_unbounded(z,M,d,alpha1,alpha2,ph,gamma,m2,m1,tmax,u0,options,alt,s);
            
            meanprey(mm) = mean(totalprey_op); minprey(mm) = min(totalprey_op); maxprey(mm) = max(totalprey_op);
            meanpred(mm) = mean(v_op(:,end)); minpred(mm) = min(v_op(:,end)); maxpred(mm) = max(v_op(:,end));
            meanmeantrait(mm) = mean(meanc_op); minmeantrait(mm) = min(meanc_op); maxmeantrait(mm) = max(meanc_op);
            meanvartrait(mm) = mean(varc_op); minvartrait(mm) = min(varc_op); maxvartrait(mm) = max(varc_op);
            wavelength(mm) = L;
            lag_prey_pred(mm) = phaselag_prey_pred; lag_pred_trait(mm) = phaselag_pred_trait; lag_mean_iqr(mm) = phaselag_mean_iqr;
        
        end
        
        save(filename,"bound_deviation","meanprey","maxprey","minprey","meanpred","minpred","maxpred","meanmeantrait","minmeantrait","maxmeantrait",...
            "meanvartrait","minvartrait","maxvartrait","wavelength","lag_pred_trait","lag_prey_pred","lag_mean_iqr")
    else
        load(filename)
    end
    
    f = figure;
    plot(bound_deviation,[meanprey;maxprey;minprey;meanpred; maxpred;minpred;meanmeantrait;maxmeantrait;minmeantrait;...
        meanvartrait;maxvartrait;minvartrait],'--o')
    grid on
    legend("Prey mean", "Prey max", "Prey min", "Pred mean", "Pred max", "Pred min", "Trait mean", "Trait max", "Trait min", "Variance mean", "Variance max", "Variance min", "Location","eastoutside")
    xlabel("Increase of trait domain size from $c \in [0,\alpha_1^{-1}]$", "Interpreter","latex")
    ylabel("Solution measures", "Interpreter","latex")
    titlestr = "";
    if m2 ~= 0.2
        titlestr = titlestr+"$, m_2="+num2str(m2)+"$";
    end
    if alpha1 ~= 0.75
        titlestr = titlestr+", $\alpha_1="+num2str(alpha1)+"$";
    end
    if alpha2 ~= 0.5
        titlestr = titlestr+", $\alpha_2="+num2str(alpha2)+"$";
    end
    if ph ~= 0.5
        titlestr = titlestr+", $p="+num2str(ph)+"$";
    end
    if gamma ~= 4
        titlestr = titlestr+", $\gamma="+num2str(gamma)+"$";
    end
    if d ~= 1e-3
        titlestr = titlestr+", $d="+num2str(d)+"$";
    end
    if titlestr ~= ""
        titlestr = eraseBetween(titlestr,1,2);
        title(titlestr,"Interpreter","latex")
    end
else

    % Get a list of all files in the current directory
    files = dir('num_sim_data/trait_bound_test_*');
    if height(files) > ceil(sqrt(height(files)))*(ceil(sqrt(height(files)))-1)
        a = ceil(sqrt(height(files))); b = a;
    else
        a = ceil(sqrt(height(files))); b = (ceil(sqrt(height(files)))-1);
    end
    f = figure;
    % Loop through each file
    cd num_sim_data\
    for i = 1:length(files)
        % Extract the filename without extension
        [~, name, ~] = fileparts(files(i).name);
        
        % Parse the filename and extract parameter values
        tokens = regexp(name, 'd([0-9dot]+)_ph([0-9dot]+)_gamma([0-9]+)_alpha1([0-9dot]+)_alpha2([0-9dot]+)_m1([0-9dot]+)_m2([0-9dot]+)', 'tokens');
        
        % Check if tokens are found
        if ~isempty(tokens)
            params = tokens{1}; % Access the matched tokens
            
            % Assign values to variables
            d = str2double(strrep(params{1}, 'dot', '.'));
            ph = str2double(strrep(params{2}, 'dot', '.'));
            gamma = str2double(params{3});
            alpha1 = str2double(strrep(params{4}, 'dot', '.'));
            alpha2 = str2double(strrep(params{5}, 'dot', '.'));
            m1 = str2double(strrep(params{6}, 'dot', '.'));
            m2 = str2double(strrep(params{7}, 'dot', '.'));
            
            % Display the variables
            fprintf('File: %s\n', files(i).name);
            fprintf('d = %.3f, ph = %.3f, gamma = %d, alpha1 = %.2f, alpha2 = %.2f, m1 = %.2f, m2 = %.2f\n', ...
                d, ph, gamma, alpha1, alpha2, m1, m2);
            
            % Load the .mat file (optional, if needed)
            % data = load(files(i).name);
        else
            fprintf('Filename %s does not match the expected format.\n', files(i).name);
        end
        load(files(i).name)
        subplot(a,b,i)

        
        % Generate 12 distinct colors using the HSV colormap
        colors = hsv(12);
        
        % Define 12 distinct marker symbols
        markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*', '+', 'x'};
        
        % Combine the 12 curves into a matrix rows. Each row represents one curve.
        data = [meanprey; maxprey; minprey; ...
                meanpred; maxpred; minpred; ...
                meanmeantrait; maxmeantrait; minmeantrait; ...
                meanvartrait; maxvartrait; minvartrait];
        
        
         hold on;
        
        % Loop through the 12 curves and plot each with its unique color and marker
        for k = 1:12
            plot(bound_deviation, data(k,:), '--', ...
                 'Color', colors(k,:), ...
                 'Marker', markers{k}, ...
                 'MarkerSize', 4, ...
                 'LineWidth', 1);
        end

        grid on
%         legend("Prey mean", "Prey max", "Prey min", "Pred mean", "Pred max", "Pred min", "Trait mean", "Trait max", "Trait min", "Variance mean", "Variance max", "Variance min", "Location","eastoutside")
%         xlabel("Increase of trait domain size from $c \in [0,\alpha_1^{-1}]$", "Interpreter","latex")
%         ylabel("Solution measures", "Interpreter","latex")
        titlestr = "";
        if m2 ~= 0.2
            titlestr = titlestr+", $m_2="+num2str(m2)+"$";
        end
        if alpha1 ~= 0.75
            titlestr = titlestr+", $\alpha_1="+num2str(alpha1)+"$";
        end
        if alpha2 ~= 0.5
            titlestr = titlestr+", $\alpha_2="+num2str(alpha2)+"$";
        end
        if ph ~= 0.5
            titlestr = titlestr+", $p="+num2str(ph)+"$";
        end
        if gamma ~= 4
            titlestr = titlestr+", $\gamma="+num2str(gamma)+"$";
        end
        if d ~= 1e-3
            titlestr = titlestr+", $d="+num2str(d)+"$";
        end
        if titlestr ~= ""
            titlestr = eraseBetween(titlestr,1,2);
            title(titlestr,"Interpreter","latex")
        end
    end
    xlabel("Increase of trait domain size from $c \in [0,\alpha_1^{-1}]$", "Interpreter","latex", "Position", [7.5,-0.1])
    ylabel("Solution measures", "Interpreter","latex", "Position", [-1,2])
    legend("Prey mean", "Prey max", "Prey min", ...
       "Pred mean", "Pred max", "Pred min", ...
       "Trait mean", "Trait max", "Trait min", ...
       "Variance mean", "Variance max", "Variance min", ...
       'Location', 'eastoutside', 'NumColumns', 2);

end

set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[10 -5 21 30])