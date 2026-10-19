%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   A Matlab program to analyze the effect of increasing/decreasing the   %
%   the persistence in the stochastic element of a simple stochastic      %
%   growth model.                                                         %
%	The stochastic element in question is a 2-state Markov process        %
%   defined at the bottom, and includes some error prevention steps.      %
%   The model itself is a simple stochastic growth model with a constant  %
%   borrowing constraint in a closed economy with no government.          %
%   By: Tomás Roda Faria Rodrigues                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('An analysis of productivity shock persistance in simple ');
disp('neoclassical growth models with borrowing constraints.');

%%%  Figure selection  %%%
% To allow for easy comparison between scenarios, you can toggle each
% figure by setting a 1 for the figures you wish to see.
showBellman = 0;    % View the Bellman equations
showPolicy = 0;     % View the policy functions
showCapDist = 1;    % View the capital distribution
showCons = 0;       % View the consumption, investment, and productivity shocks
showCapStock = 0;   % View the capital stock

%%%  Set parameter values  %%%
alpha  = 0.33;                   % Capital's share of income
beta   = 0.95;                   % Subjective discount factor 
prob1   = [ .51 .49; .49 .51 ];  % Default matrix for TFP (high-low states)
prob2   = [ .60 .40; .40 .60 ];  % More persistence in TFP
prob3   = [ .40 .60; .60 .40 ];  % More volatility in TFP
prob = { prob1 prob2 prob3 };
labels = { 'Default', 'More persistance', 'More volatility'};
delta  = .04;                    % Capital depreciation rate
A_high = 1.1;                    % "High" value for TFP
A_low  = 0.9;                    % "Low" value for TFP

%%%  Construct capital grid  %%%
mink =  4;                          % Minimum value of the capital grid
maxk =  10.0;                       % Maximum value of the capital grid   
ink  =  0.01;                       % Increment size for capital grid
nk   = round((maxk-mink)/ink+1);    % Number of grid points
kgrid = zeros(nk,0);                
for i=1:nk
   kgrid(i)=(i-1)*ink + mink; 
end

%%%  Single period return function  %%%
R_high=-10000*ones(nk,nk);  % Utility when A = high
R_low=-10000*ones(nk,nk);   % Utility when A = low

%%%  Loop through all possible current states  %%%
for i=1:nk
    k=kgrid(i);
    
    %  Loop through all possible future states
    for j=1:nk
        kprime    = kgrid(j);
        invest    = kprime - (1-delta)*k; 
        cons_high = A_high*k^(alpha) - invest;
        cons_low  = A_low*k^(alpha) - invest;
        if cons_high > 0
           R_high(j,i)= log(cons_high);
        end
        if cons_low > 0
           R_low(j,i)= log(cons_low);
        end

        %  Implement borrowing constraint
        if kprime < 6.2
           R_high(j,i)= -10000;
           R_low(j,i)= -10000;
        end
    end
end

%%%  Preparation of variables for VFI below  %%%
v       = zeros(nk,2);
decis   = zeros(nk,2);
iter    = 0;
metric1 = 10;
metric2 = 10;
[rs,cs] = size(R_high);

for g=1:numel(prob)
    p=prob{g};
    
    fprintf('\n # Scenario: %s \n', labels{g});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Iterating on Bellman's equation  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    while ((metric1 ~= 0) || (metric2 > 0.01))
        
        % For high & low productivity, calculate the best choice (k') given each
        % state variable (k) and its row location (tdecis).
        [tv_high,tdecis_high] = max(R_high + beta*p(1,1)*v(:,1)*ones(1,nk) ...
            + beta*p(1,2)*v(:,2)*ones(1,nk));
        [tv_low,tdecis_low] = max(R_low + beta*p(2,1)*v(:,1)*ones(1,nk)...
            + beta*p(2,2)*v(:,2)*ones(1,nk));
        tdecis=[tdecis_high' tdecis_low'];
        tv=[tv_high' tv_low'];
    
        % Conditions ensuring convergence:
        metric1 = max(any(tdecis-decis));       % Stop if choice = state
        metric2 = max(abs((v(:)-tv(:))./v(:))); % Stop if value functions are sufficiently close
    
        v=tv;
        decis=tdecis;
        iter=iter+1;
        % If you would like to display VFI steps, uncomment below:
        % disp([ iter metric1 metric2 ]);
    end
    
    % Count where constraint was reached
    decis=(decis-1)*ink + mink; % Turn row index into corresponding capital unit
    constraint_binds = zeros(601, 2);
    for i=1:601
        if decis(i,1) <= 6.2 + ink
            constraint_binds(i,1) = 1;
        end
        if decis(i,2) <= 6.2 + ink
            constraint_binds(i,2) = 1;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%       Model Visualizations        %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Figure 1: Value function
    if true(showBellman)
        figure(1+5*(g-1))
        plot(kgrid',v(:,1),'-',kgrid',v(:,2),':');
        title(sprintf('Stochastic Growth Model: Bellman (%s)', labels{g}));
        axis([ 4 10 0 12 ]);
        xlabel('k');
        ylabel('V(A,kt+1)');
    end
    
    % Figure 2: Policy functions (Low prod., high prod., and 45 degree line)
    if true(showPolicy)
        figure(2+5*(g-1))
        plot(kgrid',decis(:,1),'.',kgrid',decis(:,2),':',kgrid',kgrid','-', 'LineWidth',1.5);
        yline(6.2, '--r', 'k minimum', 'LineWidth',1);
        hold on;
        patch([mink maxk maxk mink], [mink mink 6.2 6.2], [0.8 0.8 0.8], ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none');
        title(sprintf('Stochastic Growth Model: Policy Functions (%s)', labels{g}));
        hold off;
        axis([ mink maxk mink maxk ]);
        xlabel('k');
        ylabel('kt+1');
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%         Transition Matrix         %%%
    %%%                 &                 %%% 
    %%%        Capital Distribution       %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Joining optimal policy choice and the 2-state Markov process
    % to create a transition matrix between periods
    g_low=zeros(cs,cs);
    g_high=zeros(cs,cs);
    for i=1:cs
        g_high(i,tdecis_high(i))=1; % For each k (column), set 1 if k' (row) is optimal
        g_low(i,tdecis_low(i))=1;
    end
   
    % Combine optimal policy matrices (g_high/g_low) with transition
    % probabilities
    trans=[ p(1,1)*g_high p(1,2)*g_high; p(2,1)*g_low p(2,2)*g_low];
    trans= trans';
    probst = (1/(2*nk))*ones(2*nk,1);
    test = 1;
    
    % Iterate on capital distribution until it converges to the steady state
    % distribution
    while test > 10^(-8)
       probst1 = trans*probst;
       test=max(abs(probst1-probst));
       probst = probst1;
    end
    
    % Calculate capital level statistics
    kk=decis(:);
    meanK=probst'*kk;                   % Mean of capital level
    varK = probst' *(kk-meanK).^2;      % Variance of capital level
    sdK = sqrt(varK);                   % Standard deviation of capital level
    
    lambda=zeros(cs,2);
    lambda(:)=probst; % Probability of having a certain (A, k) combination
    probk=sum(lambda');     
    probk=probk'; % Probability of getting capital k, regardless of shock
    
    constraint_high = (tdecis_high == 221);
    constraint_low = (tdecis_low == 221);
    bind_prob = sum(probst(1:601).*constraint_high') + ... % Probability of hitting the constraint
        sum(probst(602:1202).*constraint_low');
    
    
    disp('Capital Distribution Statistics: ');
    disp('');
    disp(' => Mean of K ');
    disp( meanK );
    disp('');
    disp(' => Standard deviation of K');
    disp( sdK )
    disp('');
    disp(' => Probability of constraint binding');
    disp( bind_prob )
    
    
    % Figure 3: Capital distribution
    if true(showCapDist)
        figure(3+5*(g-1))
        plot(kgrid,probk);
        title(sprintf('Distribution of capital: (%s)', labels{g}));
        xlabel('Capital');
        xline(6.2, '--r', 'k minimum', 'LineWidth', 1.5);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%      Life History Simulation      %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp("Simulation of agent's life");
    kmark = 200;
    k = kgrid(kmark);          % Initial level of capital/assets
    n = 500;                   % Total number of periods to simulate
    s0 = 1;                    % Initial Markov state 
    states   = zeros(n-1,2);
    controls = zeros(n-1,2);
    [chain,state] = markov(p,n,s0);
    for i = 1:n-1
        if chain(i) == 1       % High productivity shock
           kprime = decis(kmark,1);
           invest = kprime - (1-delta)*k; 
           cons   = A_high*k^(alpha) - invest; 
           kmark = tdecis(kmark,1);
        elseif chain(i) == 2   % Low productivity shock
           kprime = decis(kmark,2);
           invest = kprime - (1-delta)*k; 
           cons   = A_low*k^(alpha) - invest; 
           kmark = tdecis(kmark,2);
        else
          disp('ERROR: There is an error with the simulation chain.');
        end
        states(i,:) = [ k chain(i) ];
        controls(i,:) = [ cons kprime ];
        investment(i) = invest;
        k = kprime;
        Ktest(i) = (k == 6.2);
    end
    
    Ashocks=zeros(n-1,1);
    for i=1:n-1
        if states(i,2)==1
            Ashocks(i)=A_high;
        else
            Ashocks(i)=A_low;
        end
    end
    
    % Figure 4: Consumption, investment, and productivity shocks over time
    if true(showCons)
        figure(4+5*(g-1))
        plot((1:n-1)',investment','-',(1:n-1)',controls(:,1),(1:n-1)',Ashocks(1:n-1),':');
        title(sprintf("Agent's consumption, investments, and productivity shocks (%s)", labels{g}));
        xlabel('Period');
        ylabel('Ct, It, At');
    end
    
    % Figure 5: Capital state space and constraint binds
    constrained_periods = find((states(:,1) <= 6.2 + ink));
    
    if true(showCapStock)
        figure(5+5*(g-1))
        plot((1:n-1)',controls(:,2),'-');
        hold on
        scatter(constrained_periods, states(constrained_periods, 1), 10, 'w', 'filled');
        yline(6.2, '--r', 'LineWidth', 1.5);
        hold off
        title(sprintf("Agent's capital stock over time (%s)", labels{g}));
        axis([ 0 500 5 8 ]);
        xlabel('Period');
        ylabel('kt+1');
    end
    
    %%%  Fraction of periods where constraint was bound
    frac = sum(Ktest)/500;
    disp('Fraction of periods where constraint was bound:');
    disp('');
    disp(' => Fraction: ');
    disp( frac );

    metric1 = 10;
    metric2 = 10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      Markov Process Function                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [chain,state]=markov(T,n,s0)
%  This function runs a markov simulation, returning outputs:
%  1. chain => A simulation from a Markov chain of dimension T
%  2. state => A matrix recording the number of the realized state 
%              at time t
%
%  Included variables:
%    - T => The transition matrix
%    - n => The number of periods to simulate
%    - s0 => The initial state
%    - V => The value corresponding to each state

    [r, c]=size(T);
    V=1:r;

    %%%%%%%%%%%%%%%%%%   Error messaging   %%%%%%%%%%%%%%%%%%
    %%%  Filling in defaults for missing values  %%%
    if nargin == 1
        s0=1;
        n=100;
    end
    if nargin == 2
        s0=1;
    end
    
    %%%  Ensuring state-structure holds across periods  %%%
    if r ~= c
        disp('Error in markov function:');
        disp('   => Transition matrix must be square');
        return;
    end

    %%%  Ensuring all possible futures add up to one  %%%
    for k=1:r
        if sum(T(k,:)) ~= 1
            disp('Error in markov function')
            disp(['   => Row ',num2str(k),' does not sum up to one.']);
            disp(['   Instead, it sums to :', sum(T(k,:))]); 
            disp(['   Normalizing row ',num2str(k),'']);
            T(k,:)=T(k,:)/sum(T(k,:));
        end
    end

    %%%  Ensuring that the initial state is within bounds  %%%
    if s0 < 1 || s0 > r
        disp('Error in markov function')
        disp(['   => Initial state ',num2str(s0),' is out of range.']);
        disp('   Initial state defaulting to 1');
        s0=1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%  Markov simulation  %%%%%%%%%%%%%%%%%%
    X=rand(n-1,1);
    s=zeros(r,1);
    s(s0)=1;
    cum=T*triu(ones(size(T)));

    for k=1:length(X)
        state(:,k)=s;
        ppi=[0 s'*cum];
        s=((X(k)<=ppi(2:r+1)).*(X(k)>ppi(1:r)))';
    end
    chain=V*state;
end