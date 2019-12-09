%%% Estimate of the probability of error and the average error for a fixed value of q by using n number of (noisy) samples through runs=500 independent runs
%Requires: Synthetic Samples noisy, the number of samples n, the valueof cross-over probability q, the graph of the original Tree G, the
%original correlation matrix and the maximum allowed error gamma
%Returns: The estimated error probability and the average error of the estimated marginal distributions
function [store_mean_ssTV,store_prob] = Monte_Carlo_Iterations(noisy,n,q,runs,p,G,corr_matrix_true,gamma)
        store_mean_ssTV=zeros(1,runs);
        store_prob=zeros(1,runs);
        parfor iter=1:runs % Monte-Carlo averaging: Broadcast variables can be fixed below by predefing them, however that requires further memory. There exist a memory-running time trade-off.
            %%%Estimate the correlation matrix of noisy data
            Corr_matrix_estimate_noisy=noisy(n*(iter-1)+1:n*(iter),:)'*noisy(n*(iter-1)+1:n*(iter),:)/n;
            Corr_matrix_estimate_noisy(Corr_matrix_estimate_noisy < -(1-2*q)^2) = -(1-2*q)^2; %Apply the thresholding on the estimated correlations
            Corr_matrix_estimate_noisy(Corr_matrix_estimate_noisy > +(1-2*q)^2) = +(1-2*q)^2;

            
            %%% Find the estimated structure %%%
            [Tree_est_noisy,Cost2] = UndirectedMaximumSpanningTree (Corr_matrix_estimate_noisy);
            [rowest,colest] = find(Tree_est_noisy);
            
            %To calculate the error (ssTV) we have to find the maximum total variation between the true and the estimated pairwise marginal among all the paths of the tree
            %The shortest path algorithm finds the path of each of the p choose 2 pairs of nodes. 
            Gest = digraph(rowest,colest);
            sstv=0;
            for s=1:p
                for t=s+1:p
                    Path = shortestpath(G,s,t);
                    Pathest = shortestpath(Gest,s,t);
                    corr_prod=1;
                    corr_prod_est=1/(1-2*q)^(2*length(Pathest)-2);
                    for plc=1:length(Path)-1
                         corr_prod = corr_prod * corr_matrix_true(Path(plc),Path(plc+1));                                  
                    end    
                    for plc=1:length(Pathest)-1                
                         corr_prod_est = corr_prod_est * Corr_matrix_estimate_noisy(Pathest(plc),Pathest(plc+1));                 
                    end  
                    sstv=max(sstv,0.5*abs(corr_prod-corr_prod_est)); %By the definition of the ssTV metric
                end
            end
            store_mean_ssTV(iter)=sstv/runs; %Update the table of the average errors
            if sstv>=gamma
                store_prob(iter)=1/runs;% Update the table of the error probability, probability of the event: ssTV>=gamma
            end     
        end
end

