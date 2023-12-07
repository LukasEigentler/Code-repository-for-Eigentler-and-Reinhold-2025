function  [t,v,totalprey,medianc,meanc,L,v_op,totalprey_op,t_op,medianc_op,meanc_op,interq_trait_op,phaselag_prey_pred,phaselag_pred_trait] = prey_defence_single_run_fun(c,M,d,alpha1,alpha2,ph,gamma,m2,m1,tmax,u0,options)
% This script is a function that simulates the model once and calculates a
% number of output quantities of interest.
rerun = 1;
tvec = [0,tmax];
while rerun == 1
    [t,v] = ode15s(@(t,v) pred_prey_prey_defence_ode(v,c,M,d,alpha1,alpha2,ph,gamma,m2,m1), tvec, u0, options);
    
    %% Postprocessing
    
    %% total over time
    if M>1
        dc = c(2)-c(1);
    else
        dc = 1;
    end
    for tt = 1:length(t)
        totalprey(tt) = sum(v(tt,1:M))*dc;
        [~,tempind] = min(abs(cumsum(v(tt,1:M))/sum(v(tt,1:M))-0.5));
        medianc(tt) = c(tempind);
        meanc(tt) = sum(v(tt,1:M).*c)/sum(v(tt,1:M));
    end
    
    
    %% isolate solution over one period
    [~,tind] = min(abs(t-(tmax-100))); % find index 100 time units before end
    locmaxind = find(islocalmax(totalprey)); % find locmax of total prey
    if length(locmaxind)>3 && (max(totalprey(tind:end)) - min(totalprey(tind:end)))/mean(totalprey(tind:end)) > 1e-3 % check if there are oscillations and store info
        cycles = 1;
        

    else
        cycles = 0;
    end
    
    if cycles == 1
        startind = locmaxind(end-2); endind = locmaxind(end-1); % select second to last period
    else 
        endind = length(totalprey); startind = endind-100; % for constant sol select just a few indices so that code below works
    end
        v_op = v(startind:endind,:); % solution across one period
        totalprey_op = totalprey(startind:endind); % total prey across one period
        totalpred_op = v_op(:,end); % predator across one period
        
    %% isolate solution over the previous period
    
    if cycles == 1
        startind1 = locmaxind(end-3); endind1 = locmaxind(end-2); % select second to last period
    else 
        endind1 = length(totalprey); startind1 = endind1-100; % for constant sol select just a few indices so that code below works
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
    if cycles == 1
        L = t(endind) - t(startind); % wavelength
    else
        L = NaN; % no wavelength if no cycles
    end
    
    %% quantify solution across one period - trait interquartile range
    np = length(c);
    interq_trait_op=zeros(1,length(t_op)); traitdist = zeros(np,length(t_op));
    for i=1:length(t_op) % loop through all times
        traitdist(1:np,i) = v_op(i,1:np)/sum(v_op(i,1:np)); % frequency of trait value
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
        interq_trait_op(i) = cp75-cp25; % interquartile range of traits
    end
    
       
    
    %% quantify solution across one period - phaseshifts
    if cycles == 1 
        [~,maxpreyind] = max(totalprey_op); [~,maxpredind] = max(totalpred_op); [~,maxtraitind] = max(meanc_op); % find lcoation of max of all densities of interest
        phaselag_prey_pred = (t_op(maxpredind) - t_op(maxpreyind))/L; % phaselag between prey and predator rel to period
        phaselag_pred_trait = (t_op(maxtraitind) - t_op(maxpredind))/L; % pahselag between predator and mean trait rel to period
        % correct negative phaselags to number between 0 and 1
        if phaselag_prey_pred < 0 
            phaselag_prey_pred = 1+phaselag_prey_pred;
        end
        if phaselag_pred_trait < 0 
            phaselag_pred_trait = 1+phaselag_pred_trait;
        end
    else
        phaselag_pred_trait = NaN; % no phaseshifts if no cycles
        phaselag_prey_pred = NaN;
    end




end