function [gbestval, gbest, curve_val, curve_fes] = alg_PSO(fhd, func_num, D, NP, MaxFEs, MaxGen, lb, ub)
% ALG_PSO  PSO goc voi trong so quan tinh giam tuyen tinh (LDIW-PSO)
% Y. Shi and R. C. Eberhart, "A modified particle swarm optimizer",
% IEEE ICEC 1998.
%   w giam tuyen tinh 0.9 -> 0.4, c1 = c2 = 2, Vmax = 0.2*(ub - lb)
%
% Dau ra:
%   gbestval  - gia tri tot nhat tim duoc
%   gbest     - nghiem tot nhat (1 x D)
%   curve_val - gia tri tot nhat sau moi the he
%   curve_fes - so lan danh gia ham tuong ung voi curve_val

    w_max = 0.9; w_min = 0.4;
    c1 = 2.0; c2 = 2.0;

    lb = repmat(lb, 1, D / numel(lb));
    ub = repmat(ub, 1, D / numel(ub));
    LB = repmat(lb, NP, 1);
    UB = repmat(ub, NP, 1);
    Vmax = repmat(0.2 * (ub - lb), NP, 1);

    % Khoi tao
    pos = LB + (UB - LB) .* rand(NP, D);
    vel = -Vmax + 2 * Vmax .* rand(NP, D);
    [fit, FEs] = eval_pop(fhd, pos, func_num, 0, MaxFEs);

    pbest = pos;
    pbestval = fit;
    [gbestval, id] = min(pbestval);
    gbest = pbest(id, :);

    curve_val = zeros(1, MaxGen);
    curve_fes = zeros(1, MaxGen);
    gen = 1;
    curve_val(gen) = gbestval;
    curve_fes(gen) = FEs;

    while FEs < MaxFEs && gen < MaxGen
        gen = gen + 1;
        w = w_max - (w_max - w_min) * (gen - 1) / (MaxGen - 1);

        vel = w * vel ...
            + c1 * rand(NP, D) .* (pbest - pos) ...
            + c2 * rand(NP, D) .* (repmat(gbest, NP, 1) - pos);
        vel = max(-Vmax, min(Vmax, vel));
        pos = pos + vel;
        pos = max(LB, min(UB, pos));

        [fit, FEs] = eval_pop(fhd, pos, func_num, FEs, MaxFEs);

        improved = fit < pbestval;
        pbest(improved, :) = pos(improved, :);
        pbestval(improved) = fit(improved);
        [val, id] = min(pbestval);
        if val < gbestval
            gbestval = val;
            gbest = pbest(id, :);
        end

        curve_val(gen) = gbestval;
        curve_fes(gen) = FEs;
    end

    curve_val = curve_val(1:gen);
    curve_fes = curve_fes(1:gen);
end
