function  [t,v,totalprey,medianc,meanc,L,v_op,totalprey_op,t_op,medianc_op,meanc_op,varc_op,phaselag_prey_pred,phaselag_pred_trait,phaselag_mean_var,varc] = prey_defence_single_run_fun_unbounded(z,M,d,alpha1,alpha2,ph,gamma,m2,m1,tmax,u0,options,alt,s)
% This script is a function that simulates the model once and calculates a
% number of output quantities of interest.
rerun = 1;
tvec = [0,tmax];
% c = 1./alpha1./(1+exp(-z));
c = z;
% c(z<0) = 0; c(z>1/alpha1) = 1/alpha1;
while rerun == 1
    [t,v] = ode15s(@(t,v) pred_prey_prey_defence_ode_unbounded(v,z,M,d,alpha1,alpha2,ph,gamma,m2,m1,alt,s), tvec, u0, options);
    
    %% Postprocessing
    
    %% total over time
    if M>1
        dz = z(2)-z(1);
    else
        dz = 1;
    end
    totalprey = []; medianc = []; meanc = []; varc = [];
    for tt = 1:length(t)
        totalprey(tt) = sum(v(tt,1:M))*dz;
        [~,tempind] = min(abs(cumsum(v(tt,1:M))/sum(v(tt,1:M))-0.5));
        medianc(tt) = c(tempind);
        meanc(tt) = sum(v(tt,1:M).*c)/sum(v(tt,1:M));
        varc(tt) = sum(v(tt,1:M).*(c-meanc(tt)).^2)/sum(v(tt,1:M));
    end
    
    
    %% isolate solution over one period
    [~,tind] = min(abs(t-(tmax-100))); % find index 100 time units before end
    locminind = find(islocalmin(v(:,end))); % find locmax of total prey
    if length(locminind)>3 && (max(totalprey(tind:end)) - min(totalprey(tind:end)))/mean(totalprey(tind:end)) > 1e-3 % check if there are oscillations and store info
        cycles = 1;
        

    else
        cycles = 0;
    end
    
    if cycles == 1
        startind = locminind(end-2); endind = locminind(end-1); % select second to last period
    else 
        endind = length(totalprey); startind = endind-1; % for constant sol select just a few indices so that code below works
    end
        v_op = v(startind:endind,:); % solution across one period
        totalprey_op = totalprey(startind:endind); % total prey across one period
        totalpred_op = v_op(:,end); % predator across one period
        
    %% isolate solution over the previous period
    
    if cycles == 1
        startind1 = locminind(end-3); endind1 = locminind(end-2); % select second to last period
    else 
        endind1 = length(totalprey); startind1 = endind1-1; % for constant sol select just a few indices so that code below works
    end
    v_op1 = v(startind1:endind1,:); % solution across one period    
    totalprey_op1 = totalprey(startind1:endind1); % total prey across one period
    totalpred_op1 = v_op1(:,end); % predator across one period
    %% check if solution has reached equilibrium - if not run for longer
    if abs(max(totalpred_op1)-min(totalpred_op1) - (max(totalpred_op) - min(totalpred_op)))/(max(totalpred_op) - min(totalpred_op)) >1e-2 || ...
            abs(max(totalprey_op1)-min(totalprey_op1) - (max(totalprey_op) - min(totalprey_op)))/(max(totalprey_op) - min(totalprey_op)) >1e-2
        rerun = 1;
        u0 = v(end,:);
        tvec = [tvec(end),tvec(end)+tmax];
        fprintf("Additional run needed. Now up to "+num2str(tvec(end)) + "\n");
    else
        rerun = 0;
    end
end

    %% quantify solution across one period - densities, mean traits, wavelength
    
    
    t_op = t(startind:endind) - t(startind); % time vector of that period
    medianc_op = medianc(startind:endind); % median trait across one period
    meanc_op = meanc(startind:endind); % mean trait across one period
    varc_op = varc(startind:endind); % trait variance across one period
    if cycles == 1
        L = t(endind) - t(startind); % wavelength
    else
        L = NaN; % no wavelength if no cycles
    end
    
    %% quantify solution across one period - trait interquartile range
%     np = length(c);
%     interq_trait_op=zeros(1,length(t_op)); traitdist = zeros(np,length(t_op));
%     for i=1:length(t_op) % loop through all times
%         traitdist(1:np,i) = v_op(i,1:np)/sum(v_op(i,1:np)); % frequency of trait value
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
%         interq_trait_op(i) = sqrt(varc_op(i)); % interquartile range of traits
%     end
    

    
    %% quantify solution across one period - phaseshifts
    if cycles == 1 
        [~,maxpreyind] = max(totalprey_op); [~,maxpredind] = max(totalpred_op); 
        [~,maxtraitind] = max(meanc_op); [~,maxvarind] = max(varc_op); % find lcoation of max of all densities of interest
        phaselag_prey_pred = (t_op(maxpredind) - t_op(maxpreyind))/L; % phaselag between prey and predator rel to period
        phaselag_pred_trait = (t_op(maxtraitind) - t_op(maxpredind))/L; % pahselag between predator and mean trait rel to period
        phaselag_mean_var = (t_op(maxvarind) - t_op(maxtraitind))/L; % pahselag between mean trait and iqr trait rel to period
        % correct negative phaselags to number between 0 and 1
        if phaselag_prey_pred < 0 
            phaselag_prey_pred = 1+phaselag_prey_pred;
        end
        if phaselag_pred_trait < 0 
            phaselag_pred_trait = 1+phaselag_pred_trait;
        end
        if phaselag_mean_var < 0 
            phaselag_mean_var = 1+phaselag_mean_var;
        end
    else
        phaselag_pred_trait = NaN; % no phaseshifts if no cycles
        phaselag_prey_pred = NaN;
        phaselag_mean_var = NaN;
    end




end