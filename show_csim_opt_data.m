function show_csim_opt_data(dirLoc)
if nargin < 1
    dirLoc = './';
end
figure;

dLength = 0;
while (1) 
    d = dir([dirLoc 'population_it_*.fits']);
    if dLength ~= length(d)
        pause(0.1);
        
        subplot(2,2,1);
        fid = fopen('optValHistory.txt', 'r');
        a = fscanf(fid, '%d %f');
        c = a(1:2:end);
        v = a(2:2:end);
        semilogy(c,v,'ro-');
        grid on;
%         ylim([0.99*min(v(2:end)) 1.01*max(v(2:end))]);

        pop = fitsread('initial_population.fits');
        s1 = 1;
        s2 = 2;
        subplot(2,2,3);
        plot(pop(s1,:), pop(s2,:), '.');
        hold on; 
        ax = gca;
        ax.ColorOrderIndex = 1;
        plot(mean(pop(s1,:)), mean(pop(s2,:)), 'r*');
        
        axis equal;
        
        s1 = 400;
        s2 = 600;
        subplot(2,2,4);
        plot(pop(s1,:), pop(s2,:), '.');
        hold on; 
        ax = gca;
        ax.ColorOrderIndex = 1;
        plot(mean(pop(s1,:)), mean(pop(s2,:)), 'r*');
        axis equal;
        
        for i=1:length(d)
            pop = fitsread([dirLoc d(i).name]);
            
            s1 = 1;
            s2 = 2;
            subplot(2,2,3);
            plot(pop(s1,:), pop(s2,:), '.');
            plot(mean(pop(s1,:)), mean(pop(s2,:)), 'r*');

            s1 = 400;
            s2 = 600;
            subplot(2,2,4);
            plot(pop(s1,:), pop(s2,:), '.');
            plot(mean(pop(s1,:)), mean(pop(s2,:)), 'r*');
        end
        
        subplot(2,2,2);
        m = mean(pop, 2); % array of mean of the popoulation
        n = zeros(size(pop, 2), 1);
        for i=1:size(pop,2) 
            n(i) = norm(pop(:,i) - m); % need to do it by population element to get vector norm
        end
        hist(n,30);
        title('Distribution of distance of sag population from mean');

        dLength = length(d);
    end
    pause(1);
end