
% 1 - basic SNR TE / TR
% 2 - SNR per time unit
c_case = 4; 

T2_3T = 70e-3;
T2_7T = 45e-3;

T1_3T = 800e-3;
T1_7T = 1200e-3;

TE = linspace(40e-3,160e-3,101);
TR = linspace(0,20,101);

TR_ref = 5;

if (c_case == 1) || (c_case == 5)
    f = @(TR, T1) (1 - exp(-(TR + 1e-3) / T1));% sqrt(TR_ref ./ (TR + 1e-3));
else
    f = @(TR, T1) (1 - exp(-(TR + 1e-3) / T1)) .* sqrt(TR_ref ./ (TR + 1e-3));
end

SNR_3T = 40 * 3.0 * exp(-TE' / T2_3T) * f(TR, T1_3T);

SNR_7T = 40 * 7.0 * exp(-TE' / T2_7T) * f(TR, T1_7T);


SNR_ratio = (SNR_7T ./ SNR_3T);




msf_clf;
do_contour = 0;
switch (c_case)
    case 1
        A = SNR_7T;
        sc = 150;
        t_str = 'SNR per shot (7T)';
        s_str = '%1.0f';
        cmin = 0;
        cmax = 150;
        ctick = 0:50:150;
        
    case 2
        A = SNR_7T;
        t_str = 'SNR efficiency (7T)';
        s_str = '%1.0f';
        cmin = 0;
        cmax = 150;
    case 3
        A = SNR_3T;
        sc = 150;
        t_str = 'SNR efficiency (3T)';
        s_str = '%1.0f';
    case 4
        A = SNR_ratio;
        sc = 2;
        do_contour = 1;
        t_str = 'SNR ratio (7T/3T)';
        s_str = '%1.1f';

        cmin = 0;
        cmax = 2;
        ctick = [0 1 2];
        
end

imagesc(TE([1 end]) * 1e3, TR([end 1]), A');
hold on;
if (do_contour)
    contour(TE * 1e3, fliplr(TR), A', 1.01, 'linewidth', 2, 'color', 'red');
end

caxis([cmin cmax]);
colormap parula;
h = colorbar;
set(h, 'Ticks', ctick);
xlabel('TE [ms]');
ylabel('TR [s]');
box off;
set(gca,'tickdir','out');
set(gca,'fontsize', 18);
title(t_str);

ind = [1:25:101];
set(gca,'xtick', TE(ind) * 1e3);
set(gca,'ytick', TR(ind), 'yticklabel', num2cell(TR(ind(end:-1:1))));

set(gca,'position', [0.16 0.18 0.63 0.73]);

% add b-curves
% show an interesting line
% heating: p = u * i = r * i^2,
% g propto i
% u = r * i
for c_setup = 1:2
    
    b_max = linspace(1,10,101) * 1e9;
    
    g_max = 80e-3;
    t_rf = 8e-3;
    t_epi = 15e-3;
    n_slices = 60;
    for c = 1:numel(b_max)
        
        delta = linspace(1, 100, 1001) * 1e-3;
        g_eff = g_max;
        
        while (1)
            b = (42.6e6 * 2 * pi)^2 * delta.^2 * g_eff^2 .* ( (delta + t_rf) - delta/3);
            
            delta_min = delta(find(b > b_max(c), 1, 'first'));
            
            te_min(c) = 2 * (delta_min + t_rf/2 + t_epi);
            
            energy_produced(c) = 2*delta_min * g_eff^2;
            
            p_cool = 1.5e-4 / (1e-3 * 100); % 10 ms cooloff for b = 1 ms/um2

%             p_cool = 2.3e-4 / (1e-3 * 100); % 10 ms cooloff for b = 1 ms/um2

            t_cooloff = energy_produced(c) / p_cool;
            
            if (c_setup == 1) && (t_cooloff > te_min(c))
                g_eff = g_eff / 1.01;
                continue;
            end
            
            t_shot = max(te_min(c), t_cooloff) + t_epi * 2;
            
            tr_min(c) = n_slices * t_shot;
            break;
        end
    end
    
    x = te_min;
    y = tr_min;
    
    
    plot(x * 1e3, TR(end) - y, 'k', 'color', 'white');
    f = @(x) find( abs(b_max - x) == min(abs(b_max - x)), 1);
    ind = [f(1.0e9) f(2.5e9) f(10e9)];
    plot(x(ind) * 1e3, TR(end) - y(ind), 'o', ...
        'markerfacecolor','white', ...
        'markeredgecolor', 'white', ...
        'markersize', 29);
    
    for c = 1:numel(ind)
        [~,te_ind] = min( abs(x(ind(c)) - TE));
        [~,tr_ind] = min( abs(y(ind(c)) - TR));
        
        tmp = sprintf(s_str, (A(te_ind, tr_ind)));
        
        if (numel(tmp) == 1)
            dx = 0.5;
        elseif (numel(tmp) == 2)
            dx = 2.0;
        else
            dx = 2.3;
        end
        
        dx = dx + 2.5;
        
        text(x(ind(c)) * 1e3 - dx, TR(end) - y(ind(c)), tmp, ...
            'fontsize', 16);
    end
end


