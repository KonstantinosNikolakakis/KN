%%% 3D Plots for the probability of large ssTV (greater than gamma) and the average error of the distribution estimate
% Requires: The matrices prob , mean_sstv, the maximum value and the step-size of q (qmax, qstep) and the number of samples' batches
function visual_representation(prob,mean_sstv,qstep,qmax,points)
    [x,y]=meshgrid(0.00:qstep:qmax,1:1:points);
    
    %%% Plot the probability of the error to exceed the value gamma %%%
    figure
    mesh(y,x,prob)
    xlabel('number of samples $n\times 10^3$','Interpreter','latex')
    zlabel('probability of failure $\delta$','Interpreter','latex')
    ylabel('cross-over probability $q$','Interpreter','latex')
    xlim([1 points])
    ylim([0 qmax])
    zlim([0 1])
    colorbar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Plot the average error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    mesh(y,x,mean_sstv)
    xlim([1 points])
    ylim([0 qmax])
    xlabel('number of samples $n\times 10^3$','Interpreter','latex')
    zlabel('ssTV','Interpreter','latex')
    ylabel('cross-over probability $q$','Interpreter','latex')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%Heat map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    xx = linspace(0, qmax);
    c=9.33; % The constant C varies for different values of p,gamma, for p=31 and gamma=10^(-1.5) the value c=9.33 shows that the experimental and theoretical results match  
    yy = c*(1-2*xx).^(-4); % The trend is exactly as we proved in the theorem
    plot(xx,yy,'r','linewidth',2)
    xlabel('cross-over probability $q$','Interpreter','latex')
    ylabel('number of samples $n\times 10^3$','Interpreter','latex')
    legend({'$n=c\frac{1}{(1-2q)^4}$'},'Interpreter','latex')
    hold on
    pcolor(x,y,prob) ;
    xlim([0 qmax])
    ylim([1 points])
    colorbar
    hold on
    plot(xx,yy,'r','linewidth',2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

