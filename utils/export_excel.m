function export_excel(R, xlsx_file)
% EXPORT_EXCEL  Xuat cac bang so sanh thuat toan ra file .xlsx
%   Sheet Statistics : Mean, Std, Best, Worst, Median cua sai so cuoi
%   Sheet Mean_Std   : Mean, Std va hang cua tung thuat toan tren moi ham
%   Sheet Wilcoxon   : Rank-sum (alpha = 0.05) giua thuat toan tham chieu
%                      va cac thuat toan khac, ky hieu +/=/-
%   Sheet Friedman   : Hang trung binh va kiem dinh Friedman
%   Sheet Checkpoints: Sai so trung binh tai cac moc FEs chuan CEC
%   Sheet Time       : Thoi gian chay trung binh (giay)

    cfg = R.cfg;
    algs = cfg.algs;
    nA = numel(algs);
    funcs = cfg.funcs(:);
    nF = numel(funcs);
    fnames = arrayfun(@(f) sprintf('F%d', f), funcs, 'UniformOutput', false);
    alpha = 0.05;

    if exist(xlsx_file, 'file'), delete(xlsx_file); end

    % ---------- Statistics ----------
    T = table(fnames, 'VariableNames', {'Function'});
    means = zeros(nF, nA);
    for a = 1:nA
        E = R.(algs{a}).err;            % nF x runs
        means(:, a) = mean(E, 2);
        T.([algs{a} '_Mean'])   = mean(E, 2);
        T.([algs{a} '_Std'])    = std(E, 0, 2);
        T.([algs{a} '_Best'])   = min(E, [], 2);
        T.([algs{a} '_Worst'])  = max(E, [], 2);
        T.([algs{a} '_Median']) = median(E, 2);
    end
    writetable(T, xlsx_file, 'Sheet', 'Statistics');

    % ---------- Mean_Std + hang tren tung ham ----------
    ranks = zeros(nF, nA);
    for i = 1:nF
        ranks(i, :) = tied_rank(means(i, :))';
    end
    T = table(fnames, 'VariableNames', {'Function'});
    for a = 1:nA
        E = R.(algs{a}).err;
        T.([algs{a} '_Mean']) = mean(E, 2);
        T.([algs{a} '_Std'])  = std(E, 0, 2);
        T.([algs{a} '_Rank']) = ranks(:, a);
    end
    writetable(T, xlsx_file, 'Sheet', 'Mean_Std');

    % ---------- Wilcoxon rank-sum ----------
    ref = cfg.ref_alg;
    others = setdiff(algs, {ref}, 'stable');
    T = table(fnames, 'VariableNames', {'Function'});
    summary = cell(1, numel(others));
    for k = 1:numel(others)
        o = others{k};
        pv = zeros(nF, 1);
        sg = cell(nF, 1);
        for i = 1:nF
            x = R.(ref).err(i, :);
            y = R.(o).err(i, :);
            pv(i) = rank_sum_test(x, y);
            if pv(i) < alpha && mean(x) < mean(y)
                sg{i} = '+';            % Thuat toan tham chieu tot hon
            elseif pv(i) < alpha && mean(x) > mean(y)
                sg{i} = '-';            % Thuat toan tham chieu kem hon
            else
                sg{i} = '=';
            end
        end
        T.([ref '_vs_' o '_p']) = pv;
        T.([ref '_vs_' o '_sign']) = sg;
        summary{k} = sprintf('%d/%d/%d', sum(strcmp(sg, '+')), ...
                             sum(strcmp(sg, '=')), sum(strcmp(sg, '-')));
    end
    writetable(T, xlsx_file, 'Sheet', 'Wilcoxon');
    % Bang tong ket +/=/- dat ben duoi
    S = table(strcat(ref, '_vs_', others(:)), summary(:), ...
              'VariableNames', {'Comparison', 'Win_Tie_Loss'});
    writetable(S, xlsx_file, 'Sheet', 'Wilcoxon', 'Range', sprintf('A%d', nF + 4));

    % ---------- Friedman ----------
    avg_rank = mean(ranks, 1);
    final_rank = tied_rank(avg_rank)';
    chi2 = 12 * nF / (nA * (nA + 1)) * (sum(avg_rank .^ 2) - nA * (nA + 1) ^ 2 / 4);
    p_fr = gammainc(chi2 / 2, (nA - 1) / 2, 'upper');
    T = table(algs(:), avg_rank(:), final_rank(:), ...
              'VariableNames', {'Algorithm', 'Average_Rank', 'Final_Rank'});
    T.Friedman_Chi2 = repmat(chi2, nA, 1);
    T.Friedman_p = repmat(p_fr, nA, 1);
    writetable(T, xlsx_file, 'Sheet', 'Friedman');

    % ---------- Checkpoints chuan CEC ----------
    nC = numel(cfg.cec_points);
    rows = nF * nA;
    Fcol = cell(rows, 1); Acol = cell(rows, 1); M = zeros(rows, nC);
    k = 0;
    for i = 1:nF
        for a = 1:nA
            k = k + 1;
            Fcol{k} = fnames{i};
            Acol{k} = algs{a};
            M(k, :) = mean(reshape(R.(algs{a}).cec_err(i, :, :), cfg.runs, nC), 1);
        end
    end
    T = table(Fcol, Acol, 'VariableNames', {'Function', 'Algorithm'});
    for c = 1:nC
        T.(sprintf('FEs_%g', cfg.cec_points(c) * cfg.MaxFEs)) = M(:, c);
    end
    writetable(T, xlsx_file, 'Sheet', 'Checkpoints');

    % ---------- Time ----------
    T = table(fnames, 'VariableNames', {'Function'});
    for a = 1:nA
        T.([algs{a} '_Time_s']) = mean(R.(algs{a}).time, 2);
    end
    writetable(T, xlsx_file, 'Sheet', 'Time');
end
