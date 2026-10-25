function [bestval, bestx, curve_val, curve_fes] = alg_JADE(fhd, func_num, D, NP, MaxFEs, MaxGen, lb, ub)
% ALG_JADE  DE thich nghi voi kho luu tru ngoai (JADE)
% J. Zhang and A. C. Sanderson, "JADE: Adaptive Differential Evolution With
% Optional External Archive", IEEE Trans. Evol. Comput., 13(5), 945-958, 2009.
%   Chien luoc DE/current-to-pbest/1/bin co kho luu tru
%   p = 0.05, c = 0.1, mu_CR = mu_F = 0.5 ban dau, |A| <= NP

    p = 0.05;
    c = 0.1;
    mu_CR = 0.5;
    mu_F = 0.5;

    lb = repmat(lb, 1, D / numel(lb));
    ub = repmat(ub, 1, D / numel(ub));
    LB = repmat(lb, NP, 1);
    UB = repmat(ub, NP, 1);

    pop = LB + (UB - LB) .* rand(NP, D);
    [fit, FEs] = eval_pop(fhd, pop, func_num, 0, MaxFEs);
    [bestval, id] = min(fit);
    bestx = pop(id, :);

    archive = zeros(0, D);
    n_pbest = max(1, round(p * NP));

    curve_val = zeros(1, MaxGen);
    curve_fes = zeros(1, MaxGen);
    gen = 1;
    curve_val(gen) = bestval;
    curve_fes(gen) = FEs;

    while FEs < MaxFEs && gen < MaxGen
        gen = gen + 1;

        % Sinh CR_i ~ N(mu_CR, 0.1) cat ve [0, 1]
        CR = min(max(mu_CR + 0.1 * randn(NP, 1), 0), 1);

        % Sinh F_i ~ Cauchy(mu_F, 0.1): > 1 thi cat ve 1, <= 0 thi sinh lai
        Fi = mu_F + 0.1 * tan(pi * (rand(NP, 1) - 0.5));
        bad = Fi <= 0;
        while any(bad)
            Fi(bad) = mu_F + 0.1 * tan(pi * (rand(sum(bad), 1) - 0.5));
            bad = Fi <= 0;
        end
        Fi = min(Fi, 1);

        % Chon x_pbest tu top 100p% ca the tot nhat
        [~, order] = sort(fit);
        pb = order(randi(n_pbest, NP, 1));

        % r1 tu P \ {i}, r2 tu (P U A) \ {i, r1}
        union_pop = [pop; archive];
        NU = size(union_pop, 1);
        r1 = zeros(NP, 1);
        r2 = zeros(NP, 1);
        for i = 1:NP
            r1(i) = randi(NP - 1);
            if r1(i) >= i, r1(i) = r1(i) + 1; end
            r2(i) = randi(NU);
            while r2(i) == i || r2(i) == r1(i)
                r2(i) = randi(NU);
            end
        end

        V = pop + Fi .* (pop(pb, :) - pop) + Fi .* (pop(r1, :) - union_pop(r2, :));
        V = repair_bounds(V, pop, LB, UB);

        % Lai ghep nhi thuc
        mask = rand(NP, D) < CR;
        jrand = randi(D, NP, 1);
        mask(sub2ind([NP, D], (1:NP)', jrand)) = true;
        U = pop;
        U(mask) = V(mask);

        [ufit, FEs] = eval_pop(fhd, U, func_num, FEs, MaxFEs);

        % Chon loc va cap nhat kho luu tru
        success = ufit < fit;
        archive = [archive; pop(success, :)]; %#ok<AGROW>
        if size(archive, 1) > NP
            keep = randperm(size(archive, 1), NP);
            archive = archive(keep, :);
        end
        pop(success, :) = U(success, :);
        fit(success) = ufit(success);

        % Cap nhat mu_CR (trung binh so hoc) va mu_F (trung binh Lehmer)
        if any(success)
            SCR = CR(success);
            SF = Fi(success);
            mu_CR = (1 - c) * mu_CR + c * mean(SCR);
            mu_F = (1 - c) * mu_F + c * sum(SF .^ 2) / sum(SF);
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
