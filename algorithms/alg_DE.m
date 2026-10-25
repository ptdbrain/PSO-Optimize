function [bestval, bestx, curve_val, curve_fes] = alg_DE(fhd, func_num, D, NP, MaxFEs, MaxGen, lb, ub)
% ALG_DE  Tien hoa vi phan goc DE/rand/1/bin
% R. Storn and K. Price, "Differential Evolution - A Simple and Efficient
% Heuristic for Global Optimization over Continuous Spaces",
% J. Global Optimization, 11, 341-359, 1997.
%   F = 0.5, CR = 0.9
%   Xu ly bien: dat lai ve trung diem giua ca the cha va bien bi vi pham.

    F = 0.5;
    CR = 0.9;

    lb = repmat(lb, 1, D / numel(lb));
    ub = repmat(ub, 1, D / numel(ub));
    LB = repmat(lb, NP, 1);
    UB = repmat(ub, NP, 1);

    pop = LB + (UB - LB) .* rand(NP, D);
    [fit, FEs] = eval_pop(fhd, pop, func_num, 0, MaxFEs);
    [bestval, id] = min(fit);
    bestx = pop(id, :);

    curve_val = zeros(1, MaxGen);
    curve_fes = zeros(1, MaxGen);
    gen = 1;
    curve_val(gen) = bestval;
    curve_fes(gen) = FEs;

    while FEs < MaxFEs && gen < MaxGen
        gen = gen + 1;

        % Dot bien: chon r1, r2, r3 khac nhau va khac i
        r = zeros(NP, 3);
        for i = 1:NP
            idx = randperm(NP - 1, 3);
            idx(idx >= i) = idx(idx >= i) + 1;
            r(i, :) = idx;
        end
        V = pop(r(:, 1), :) + F * (pop(r(:, 2), :) - pop(r(:, 3), :));
        V = repair_bounds(V, pop, LB, UB);

        % Lai ghep nhi thuc
        mask = rand(NP, D) < CR;
        jrand = randi(D, NP, 1);
        mask(sub2ind([NP, D], (1:NP)', jrand)) = true;
        U = pop;
        U(mask) = V(mask);

        [ufit, FEs] = eval_pop(fhd, U, func_num, FEs, MaxFEs);

        % Chon loc
        better = ufit <= fit;
        pop(better, :) = U(better, :);
        fit(better) = ufit(better);

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
