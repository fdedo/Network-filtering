%%% this function compute the p-value for a polya filter

function b = polya_filter(A, a, alpha, apr_lvl)
% A is adiacency matrix
% a is the free parameter of the polya method (howm many balls to reinsert
% in the urn) 
% alpha is the significance level
% apr_lvl: is a constant defining the regime under which the approximate form of the p-value (Eq. (6) of the paper below) can be used.
%     For example, setting apr_lvl = 10 will trigger the use of the approximate form for every link such that s > 10*k/a, w > 10, and s-w>10*k/a. 
%     The approximate form of the p-value is much faster to compute. 

    k_in = full(sum(A > 0));    % IN-Degree sequence
    k_out = full(sum(A' > 0));  % OUT-Degree sequence
    s_in = full(sum(A));        % IN-Strength sequence
    s_out = full(sum(A'));      % OUT-Strength sequence

    % Finding indices of non-zero entries in A (i.e., links)
    [ind1,ind2] = find(A > 0); 

    b = []; % Empty array to store links in backbone

    for i = 1:length(ind1) % loop on links

        w = A(ind1(i), ind2(i)); % weight of the current link

        %calculate p-values
        p_out = polya_cdf(w, s_out(ind1(i)), k_out(ind1(i)), a, apr_lvl);
        p_in = polya_cdf(w, s_in(ind2(i)), k_in(ind2(i)), a , apr_lvl);

        % Choose the p-value to check on, first take care of not having a
        % 1-degree node you're checking on
        
        if k_in(ind2(i)) == 1
            P = p_out;
        elseif k_out(ind1(i)) == 1
            P = p_in;
        else 
            P = min(p_in, p_out); %keep the minor one 
        end
        
        % If the p-value falls below the significance level in input, the
        % corresponding link is stored in the backbone
        if P < alpha
           b = [b; ind1(i) ind2(i)];
        end
       
    end

end

% function to compute the cumulative distribution (from x to infty)
function [p] = polya_cdf(w,s,k,a,L)
p = nan(length(w),1);
if a==0 %binomial case
    p = binocdf(w,s,1./k);
else
    %check where approximation can be performed
    idx1 = s-w >= L*((k-1)./a+1); % array dimensione w con 1 se cond è soddisfatta e 0 se non lo è
    idx2 = w>=L*max(1./a,1);
    idx3 = s>=L*max(k./a,1);
    idx4 = k>=L*(a-1+1e-20);
    idx = (double(idx1)+double(idx2)+double(idx3)+double(idx4))==4;

    %calculate the p-values that can be approximated
    p(idx) = 1/gamma(1/a)*(1-w(idx)./s(idx)).^((k(idx)-1)/a).*(w(idx).*k(idx)./(s(idx)*a)).^(1/a-1);
    
    %calculate the p-values that cannot be aprroximated
    idx = find(double(idx)==0); % trova gli indici con valore = 0
    if isempty(idx)==0 %controlla che esista almeno un indice con valore = 0
       
        for ii=1:length(idx)
            n = s(idx(ii));
            A = 1/a;
            B = (k(idx(ii))-1)./a;
            x = (0:1:w(idx(ii))-1)'; % tutti i valori da 0 a w 
            p(idx(ii)) = 1 - sum(exp(gammaln(n+1)+betaln(x+A,n-x+B)-gammaln(x+1)-gammaln(n-x+1)-betaln(A,B)));
        end 
    end
end
%catch rounding errors (use the mp-toolbox for higher precision)
p(p<0) = 0;
end

function [a_best,err] = get_ML_estimate(W)
%This function can be used to get the maximum likelihood estimation of the
%free parameter "a". It will set the filter on the network's own
%heterogeneity and it will therefore produce ultra-sparse backbones
%   W: is the adjecency matrix

    %check if the network is symmetric (i.e. undirected) and get the edge list for both cases
    if issymmetric(W)==1
        U = triu(W);
        [i,j,w] = find(U); 
    else
        [i,j,w] = find(W); 
    end
    
    %get the degrees and strengths
    k_in(:,1) = full(sum(W~=0,1)); %sum over W not equal to 0
    s_in(:,1) = full(sum(W,1));
    k_out(:,1) = full(sum(W~=0,2));
    s_out(:,1) = full(sum(W,2));
    
    w = [w;w];
    k = [k_out(i);k_in(j)];
    s = [s_out(i);s_in(j)];
    
    %get rid of the links with degree 1
    w(k==1) = [];
    s(k==1) = [];
    k(k==1) = [];
    
    f = @(a)eq_from_lhood([w k s],a);
    g = @(a)lhood([w k s],a);
    options_lsq = optimoptions(@lsqnonlin,'Display','off');
    
    %first find the maximum value of the likelihood
    %Solves nonlinear least-squares curve fitting problems of the form square
    %norm of a function 
    % it starts at the point x0=0.5 and finds a minimum of the sum of squares 
    % of the functions described in g.
    % the solution is always in the range lb=0 ≤ x ≤ ub=15
    x0 = lsqnonlin(g,0.5,0,15,options_lsq);
    
    %use this value to calculate the value that put the derivative to 0
    [a_best,err] = lsqnonlin(f,x0,0,15,options_lsq);
    
    %try higher precision if feval is not close to 0
    if err>1e-6
        options_lsq = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',2000,'FunctionTolerance',1e-15,'StepTolerance',1e-25,'OptimalityTolerance',1e-25,'Display','off');
        [a_best,err] = lsqnonlin(f,x0,0,15,options_lsq);
        if err>1e-6
            disp('Try stricter minimization options')
        end
    end
end


function [out] = eq_from_lhood( X, a )
%derivative of the likelihood that must be put equal to 0
    w = X(:,1);
    k = X(:,2);
    s = X(:,3);
    DL = a.^(-2).*(psi(a.^(-1))+((-1)+k).*psi(a.^(-1).*((-1)+k))+(-1).*k.*psi(a.^(-1).*k)+k.*psi(a.^(-1).*k+s)+(-1).*psi(a.^(-1)+w)+(-1).*((-1)+k).*psi(a.^(-1).*((-1)+k+a.*s+(-1).*a.*w)));
    DL(isnan(DL)) = 0;
    out = double(sum(DL));
    if isinf(out)==1
        out=1e100;
    end
end


function [out] = lhood( X, a )
%compute the likelihood to maximise. It's - log the product of the
%probability of each link based on the polya distribution 
% X is a matrix with weights, degrees, strengths in three coloums

    w = X(:, 1);
    k = X(:, 2);
    s = X(:, 3);
    
    P = polya_pdf(w, s, k, a);  %compute probability
    P(P==0)=1e-20;              %set to almost zero when 0
    out = sum(-log(P));         % compute the - loglikelihood

    if isinf(out)==1
        out = sign(out)*1e100;
    end
end


function [p] = polya_pdf(w,s,k,a)
%pdf of the distribution
    if a==0 %binomial case
        p = binocdf(w,s,1./k);
    else       
        n = s;
        A = 1/a;
        B = (k-1)./a;
        x = w;
        p = exp(gammaln(n+1)+betaln(x+A,n-x+B)-gammaln(x+1)-gammaln(n-x+1)-betaln(A,B));
    end
end


%%% DOMANDE
%%% perché non usare direttamente il comando sparse? (nel codice originale)
%%% cosa signfica k_in = k_in(j); ??
%%% perché le approssimazioni della cdf sono quelle 4? nel paper non le ho
%%% trovate
