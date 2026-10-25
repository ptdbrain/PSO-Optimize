function plot_convergence(R, fig_dir)
% PLOT_CONVERGENCE  Ve duong cong hoi tu (sai so trung binh qua cac lan chay)
% cho tung ham va mot hinh tong hop tat ca cac ham.

    cfg = R.cfg;
    algs = cfg.algs;
    nA = numel(algs);
    funcs = cfg.funcs;
    nF = numel(funcs);
    x = R.cp_fes;
    floor_val = 1e-8;           % Sai so 0 ve tren truc log tai 1e-8

    colors = [0.000 0.447 0.741; 0.850 0.325 0.098; ...
              0.466 0.674 0.188; 0.494 0.184 0.556; ...
              0.929 0.694 0.125; 0.301 0.745 0.933];
    markers = {'o', 's', '^', 'd', 'v', 'p'};
    mk_idx = round(linspace(1, numel(x), 10));

    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

    % Tung ham mot hinh
    for i = 1:nF
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 640 480]);
        hold on;
        for a = 1:nA
            y = mean(reshape(R.(algs{a}).curve(i, :, :), cfg.runs, []), 1);
            semilogy(x, max(y, floor_val), 'LineWidth', 1.6, ...
                'Color', colors(mod(a - 1, 6) + 1, :), ...
                'Marker', markers{mod(a - 1, 6) + 1}, 'MarkerIndices', mk_idx);
        end
        set(gca, 'YScale', 'log');
        grid on; box on;
        xlabel('So lan danh gia ham (FEs)');
        ylabel('Sai so trung binh  f(x) - f(x^*)');
        title(sprintf('CEC2017 F%d (D = %d, %d lan chay)', funcs(i), cfg.D, cfg.runs));
        legend(algs, 'Location', 'northeast');
        print(fig, fullfile(fig_dir, sprintf('F%02d.png', funcs(i))), '-dpng', '-r200');
        close(fig);
    end

    % Hinh tong hop
    ncol = 5;
    nrow = ceil(nF / ncol);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 2000 380 * nrow]);
    for i = 1:nF
        subplot(nrow, ncol, i);
        hold on;
        for a = 1:nA
            y = mean(reshape(R.(algs{a}).curve(i, :, :), cfg.runs, []), 1);
            semilogy(x, max(y, floor_val), 'LineWidth', 1.2, ...
                'Color', colors(mod(a - 1, 6) + 1, :));
        end
        set(gca, 'YScale', 'log', 'FontSize', 8);
        grid on; box on;
        title(sprintf('F%d', funcs(i)));
        if i == 1, legend(algs, 'Location', 'best', 'FontSize', 7); end
    end
    print(fig, fullfile(fig_dir, 'All_functions.png'), '-dpng', '-r150');
    close(fig);
end
