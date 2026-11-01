function [bestval, bestx, curve_val, curve_fes] = alg_LSHADE(fhd, func_num, D, MaxFEs, lb, ub)
% ALG_LSHADE  SHADE voi giam kich thuoc quan the tuyen tinh (L-SHADE)
% R. Tanabe and A. S. Fukunaga, "Improving the Search Performance of SHADE
% Using Linear Population Size Reduction", IEEE CEC 2014, 1658-1665.
% (Thuat toan thang cuoc thi CEC2014)
%
% Chuyen tu code C++ goc cua tac gia (lshade.cc, ban 1.0, nop cho cuoc thi
% CEC2014, kho CEC2014 cua P. N. Suganthan). Giu nguyen tung buoc:
%   - NP_init = round(18*D), NP_min = 4, H = 6, p = 0.11, ti le archive 2.6
%   - DE/current-to-pbest/1/bin co archive, xu ly bien trung diem voi cha
%   - Bo nho M_F, M_CR cap nhat bang trung binh Lehmer co trong so |delta f|
%   - Giam quan the tuyen tinh theo so FEs, loai ca the kem nhat
% Hai chi tiet cua code goc duoc giu lai:
%   - Archive luu ca the SAU khi da duoc thay bang con thanh cong
%     (pop(i) da gan bang con truoc khi chep vao archive).
%   - M_CR = -1 (CR = 0) chi khi tong trong so cua S_CR bang 0.

    NP_init = round(18 * D);
    NP_min = 4;
    H = 6;
    p_best_rate = 0.11;
    arc_rate = 2.6;

    lb = repmat(lb, 1, D / numel(lb));
    ub = repmat(ub, 1, D / numel(ub));

    NP = NP_init;
    pop = repmat(lb, NP, 1) + repmat(ub - lb, NP, 1) .* rand(NP, D);
    [fit, FEs] = eval_pop(fhd, pop, func_num, 0, MaxFEs);
    [bestval, id] = min(fit);
    bestx = pop(id, :);

    arc_size = round(NP * arc_rate);
    archive = zeros(arc_size, D);
    arc_count = 0;

    M_F = 0.5 * ones(H, 1);
    M_CR = 0.5 * ones(H, 1);
    mem_pos = 1;
    p_num = round(NP * p_best_rate);

    curve_val = zeros(1, MaxFEs);
    curve_fes = zeros(1, MaxFEs);
    gen = 1;
    curve_val(gen) = bestval;
    curve_fes(gen) = FEs;

    while FEs < MaxFEs
        gen = gen + 1;
        [~, sorted] = sort(fit);

        % Sinh F_i, CR_i tu mot o nho ngau nhien
        ri = randi(H, NP, 1);
        mu_f = M_F(ri);
        mu_cr = M_CR(ri);
        CR = min(max(mu_cr + 0.1 * randn(NP, 1), 0), 1);
        CR(mu_cr == -1) = 0;
        F = mu_f + 0.1 * tan(pi * (rand(NP, 1) - 0.5));
        bad = F <= 0;
        while any(bad)
            F(bad) = mu_f(bad) + 0.1 * tan(pi * (rand(sum(bad), 1) - 0.5));
            bad = F <= 0;
        end
        F = min(F, 1);

        % current-to-pbest/1 voi archive
        pb = sorted(randi(p_num, NP, 1));
        union_pop = [pop; archive(1:arc_count, :)];
        r1 = zeros(NP, 1);
        r2 = zeros(NP, 1);
        for i = 1:NP
            r1(i) = randi(NP);
            while r1(i) == i
                r1(i) = randi(NP);
            end
            r2(i) = randi(NP + arc_count);
            while r2(i) == i || r2(i) == r1(i)
                r2(i) = randi(NP + arc_count);
            end
        end
        V = pop + F .* (pop(pb, :) - pop) + F .* (pop(r1, :) - union_pop(r2, :));

        % Lai ghep nhi thuc, roi xu ly bien (trung diem voi cha)
        mask = rand(NP, D) < CR;
        jrand = randi(D, NP, 1);
        mask(sub2ind([NP, D], (1:NP)', jrand)) = true;
        U = pop;
        U(mask) = V(mask);
        U = repair_bounds(U, pop, repmat(lb, NP, 1), repmat(ub, NP, 1));

        [ufit, FEs] = eval_pop(fhd, U, func_num, FEs, MaxFEs);

        % Thay the the he
        equal = ufit == fit;
        success = ufit < fit;
        dif = abs(fit(success) - ufit(success));
        S_F = F(success);
        S_CR = CR(success);
        replace = equal | success;
        pop(replace, :) = U(replace, :);
        fit(replace) = ufit(replace);

        % Cap nhat archive (giu dung thu tu va cach thay the cua code goc)
        if arc_size > 1
            for i = find(success)'
                if arc_count < arc_size
                    arc_count = arc_count + 1;
                    archive(arc_count, :) = pop(i, :);
                else
                    archive(randi(arc_size), :) = pop(i, :);
                end
            end
        end

        % Cap nhat bo nho bang trung binh Lehmer co trong so
        if any(success)
            wk = dif / sum(dif);
            M_F(mem_pos) = sum(wk .* S_F .^ 2) / sum(wk .* S_F);
            if sum(wk .* S_CR) == 0
                M_CR(mem_pos) = -1;
            else
                M_CR(mem_pos) = sum(wk .* S_CR .^ 2) / sum(wk .* S_CR);
            end
            mem_pos = mod(mem_pos, H) + 1;
        end

        % Giam kich thuoc quan the tuyen tinh (LPSR)
        plan_NP = round((NP_min - NP_init) / MaxFEs * FEs + NP_init);
        if NP > plan_NP
            new_NP = max(plan_NP, NP_min);
            [~, order] = sort(fit);
            keep = sort(order(1:new_NP));   % Loai ca the kem nhat, giu thu tu
            pop = pop(keep, :);
            fit = fit(keep);
            NP = new_NP;

            arc_size = floor(NP * arc_rate);
            arc_count = min(arc_count, arc_size);
            p_num = round(NP * p_best_rate);
            if p_num <= 1, p_num = 2; end
        end

        [val, id] = min(fit);
        if val < bestval
            bestval = val;
            bestx = pop(id, :);
        end

        curve_val(gen) = bestval;
        curve_fes(gen) = FEs;
    end

    curve_val = curve_val(1:gen);
    curve_fes = curve_fes(1:gen);
end
