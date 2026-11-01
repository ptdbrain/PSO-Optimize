function [bestval, bestx, curve_val, curve_fes] = alg_INFO(fhd, func_num, D, MaxFEs, lb, ub)
% ALG_INFO  Toi uu dua tren trung binh co trong so cua cac vector (INFO)
% I. Ahmadianfar, A. A. Heidari, S. Noshadian, H. Chen, A. H. Gandomi,
% "INFO: An efficient optimization algorithm based on weighted mean of
% vectors", Expert Systems with Applications, 195, 116516, 2022.
%
%   Ba giai doan cho moi ca the:
%   1. Quy tac cap nhat (updating rule) dua tren MeanRule - trung binh co
%      trong so cua cac vector hieu.
%   2. Ket hop vector (vector combining).
%   3. Tim kiem cuc bo (local search) quanh nghiem tot nhat.
%   NP = 30 (nhu Main.m cua tac gia), cap nhat bat dong bo.
%
% Chuyen the tu code goc INFO.m cua tac gia (trang du an INFO, ban cap nhat
% 01-10-2022). Chi thay doi phan giao tiep: dem FEs, ghi duong hoi tu,
% alpha tinh theo FEs/MaxFEs thay cho it/MaxIt, unifrnd(0.1, 0.5) viet lai
% bang rand de khong can Statistics Toolbox, vong lap tung chieu o buoc ket
% hop vector viet dang vector hoa (tuong duong).

    NP = 30;

    lb = repmat(lb, 1, D / numel(lb));
    ub = repmat(ub, 1, D / numel(ub));

    X = repmat(lb, NP, 1) + repmat(ub - lb, NP, 1) .* rand(NP, D);
    [Cost, FEs] = eval_pop(fhd, X, func_num, 0, MaxFEs);

    [~, ind] = sort(Cost);
    bestx = X(ind(1), :);       bestval = Cost(ind(1));
    worst_x = X(ind(end), :);   worst_val = Cost(ind(end));
    I = randi([2 5]);
    better_x = X(ind(I), :);    better_val = Cost(ind(I));

    curve_val = zeros(1, MaxFEs);
    curve_fes = zeros(1, MaxFEs);
    gen = 1;
    curve_val(gen) = bestval;
    curve_fes(gen) = FEs;

    while FEs < MaxFEs
        gen = gen + 1;
        alpha = 2 * exp(-4 * FEs / MaxFEs);
        M_best = bestval; M_better = better_val; M_worst = worst_val;

        for i = 1:NP
            if FEs >= MaxFEs, break; end

            % --- 1. Quy tac cap nhat ---
            del = 2 * rand * alpha - alpha;
            sigm = 2 * rand * alpha - alpha;

            A1 = randperm(NP);
            A1(A1 == i) = [];
            a = A1(1); b = A1(2); c = A1(3);
            epsi = 1e-25 * rand;

            % Trung binh co trong so cua 3 vector ngau nhien
            omg = max([Cost(a) Cost(b) Cost(c)]);
            MM = [Cost(a) - Cost(b), Cost(a) - Cost(c), Cost(b) - Cost(c)];
            W = cos(MM + pi) .* exp(-MM / omg);
            WM1 = del * (W(1) * (X(a, :) - X(b, :)) + W(2) * (X(a, :) - X(c, :)) ...
                       + W(3) * (X(b, :) - X(c, :))) / (sum(W) + 1) + epsi;

            % Trung binh co trong so cua nghiem tot nhat, tot hon, kem nhat
            % (MM(2) = M_best - M_better giong code goc, bai bao ghi M_best - M_worst)
            omg = max([M_best M_better M_worst]);
            MM = [M_best - M_better, M_best - M_better, M_better - M_worst];
            W = cos(MM + pi) .* exp(-MM / omg);
            WM2 = del * (W(1) * (bestx - better_x) + W(2) * (bestx - worst_x) ...
                       + W(3) * (better_x - worst_x)) / (sum(W) + 1) + epsi;

            r = 0.1 + 0.4 * rand;
            MeanRule = r * WM1 + (1 - r) * WM2;

            if rand < 0.5
                z1 = X(i, :) + sigm * (rand * MeanRule) + randn * (bestx - X(a, :)) / (M_best - Cost(a) + 1);
                z2 = bestx + sigm * (rand * MeanRule) + randn * (X(a, :) - X(b, :)) / (Cost(a) - Cost(b) + 1);
            else
                z1 = X(a, :) + sigm * (rand * MeanRule) + randn * (X(b, :) - X(c, :)) / (Cost(b) - Cost(c) + 1);
                z2 = better_x + sigm * (rand * MeanRule) + randn * (X(a, :) - X(b, :)) / (Cost(a) - Cost(b) + 1);
            end

            % --- 2. Ket hop vector ---
            mu = 0.05 * randn(1, D);
            pick = rand(1, D);
            keep = rand(1, D) >= 0.5;
            u = z2 + mu .* abs(z1 - z2);
            use_z1 = pick < 0.5;
            u(use_z1) = z1(use_z1) + mu(use_z1) .* abs(z1(use_z1) - z2(use_z1));
            u(keep) = X(i, keep);

            % --- 3. Tim kiem cuc bo ---
            if rand < 0.5
                L = rand < 0.5;
                v1 = (1 - L) * 2 * rand + L;
                v2 = rand * L + (1 - L);
                Xavg = (X(a, :) + X(b, :) + X(c, :)) / 3;
                phi = rand;
                Xrnd = phi * Xavg + (1 - phi) * (phi * better_x + (1 - phi) * bestx);
                Randn = L * randn(1, D) + (1 - L) * randn;
                if rand < 0.5
                    u = bestx + Randn .* (MeanRule + randn * (bestx - X(a, :)));
                else
                    u = Xrnd + Randn .* (MeanRule + randn * (v1 * bestx - v2 * Xrnd));
                end
            end

            % Dua ve mien tim kiem (cat bien) va danh gia
            u = min(max(u, lb), ub);
            new_cost = feval(fhd, u', func_num);
            FEs = FEs + 1;

            if new_cost < Cost(i)
                X(i, :) = u;
                Cost(i) = new_cost;
                if new_cost < bestval
                    bestx = u;
                    bestval = new_cost;
                end
            end
        end

        % Cap nhat nghiem kem nhat va nghiem "tot hon"
        [~, ind] = sort(Cost);
        worst_x = X(ind(end), :);   worst_val = Cost(ind(end));
        I = randi([2 5]);
        better_x = X(ind(I), :);    better_val = Cost(ind(I));

        curve_val(gen) = bestval;
        curve_fes(gen) = FEs;
    end

    curve_val = curve_val(1:gen);
    curve_fes = curve_fes(1:gen);
end
