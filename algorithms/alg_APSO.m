function [gbestval, gbest, curve_val, curve_fes] = alg_APSO(fhd, func_num, D, NP, MaxFEs, MaxGen, lb, ub)
% ALG_APSO  Toi uu bay dan thich nghi (APSO) - phien ban goc
% Z.-H. Zhan, J. Zhang, Y. Li, H. S.-H. Chung, "Adaptive Particle Swarm
% Optimization", IEEE Trans. SMC-B, 39(6), 1362-1381, 2009.
%
% Cac thanh phan:
%   1. Uoc luong trang thai tien hoa (ESE): he so tien hoa f, phan loai mo
%      vao 4 trang thai S1 (tham do), S2 (khai thac), S3 (hoi tu),
%      S4 (nhay thoat) bang ham thuoc mo.
%   2. w(f) = 1 / (1 + 1.5 exp(-2.6 f)), gioi han [0.4, 0.9].
%   3. Dieu chinh c1, c2 theo trang thai, delta ~ U[0.05, 0.1],
%      c1, c2 trong [1.5, 2.5], c1 + c2 <= 4.
%   4. Hoc uu tu (ELS) tren gbest khi o trang thai hoi tu S3,
%      sigma giam tuyen tinh 1.0 -> 0.1.

    sigma_max = 1.0; sigma_min = 0.1;

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
    [gbestval, gid] = min(pbestval);
    gbest = pbest(gid, :);

    c1 = 2.0; c2 = 2.0;
    state = 1;                      % Trang thai ban dau: S1 (tham do)

    curve_val = zeros(1, MaxGen);
    curve_fes = zeros(1, MaxGen);
    gen = 1;
    curve_val(gen) = gbestval;
    curve_fes(gen) = FEs;

    while FEs < MaxFEs && gen < MaxGen
        gen = gen + 1;

        % --- 1. ESE: he so tien hoa f ---
        d = mean_distances(pos);
        [~, gid] = min(pbestval);
        dmin = min(d); dmax = max(d);
        if dmax - dmin < realmin
            f = 0;
        else
            f = (d(gid) - dmin) / (dmax - dmin);
        end

        state = classify_state(f, state);

        % --- 2. Trong so quan tinh ---
        w = 1 / (1 + 1.5 * exp(-2.6 * f));

        % --- 3. He so gia toc ---
        delta = 0.05 + 0.05 * rand();
        switch state
            case 1  % S1 Tham do: tang c1, giam c2
                c1 = c1 + delta;        c2 = c2 - delta;
            case 2  % S2 Khai thac: tang nhe c1, giam nhe c2
                c1 = c1 + 0.5 * delta;  c2 = c2 - 0.5 * delta;
            case 3  % S3 Hoi tu: tang nhe c1, tang nhe c2
                c1 = c1 + 0.5 * delta;  c2 = c2 + 0.5 * delta;
            case 4  % S4 Nhay thoat: giam c1, tang c2
                c1 = c1 - delta;        c2 = c2 + delta;
        end
        c1 = min(max(c1, 1.5), 2.5);
        c2 = min(max(c2, 1.5), 2.5);
        if c1 + c2 > 4.0
            s = c1 + c2;
            c1 = 4.0 * c1 / s;
            c2 = 4.0 * c2 / s;
        end

        % --- 4. ELS: hoc uu tu khi hoi tu ---
        if state == 3 && FEs < MaxFEs
            P = gbest;
            k = randi(D);
            sigma = sigma_max - (sigma_max - sigma_min) * FEs / MaxFEs;
            P(k) = P(k) + (ub(k) - lb(k)) * sigma * randn();
            P(k) = min(max(P(k), lb(k)), ub(k));
            vP = feval(fhd, P', func_num);
            FEs = FEs + 1;

            if vP < gbestval
                gbestval = vP;
                gbest = P;
            else
                % Thay the hat co pbest kem nhat bang P
                [~, wid] = max(pbestval);
                pos(wid, :) = P;
                pbest(wid, :) = P;
                pbestval(wid) = vP;
            end
        end

        % --- 5. Cap nhat van toc va vi tri ---
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

function d = mean_distances(pos)
% Khoang cach Euclid trung binh tu moi hat toi cac hat con lai
    NP = size(pos, 1);
    sq = sum(pos .^ 2, 2);
    dist2 = max(sq + sq' - 2 * (pos * pos'), 0);
    d = sum(sqrt(dist2), 2) / (NP - 1);
end

function state = classify_state(f, prev)
% Phan loai mo trang thai tien hoa (Zhan 2009, Hinh 4).
% Khi f thuoc hai trang thai cung luc, giu trang thai truoc neu no nam
% trong cap do (on dinh); neu khong, chon trang thai gan trang thai truoc
% nhat theo chu trinh S1 -> S2 -> S3 -> S4 -> S1.
    mu = zeros(1, 4);

    % S1 - Tham do
    if f <= 0.4,      mu(1) = 0;
    elseif f <= 0.6,  mu(1) = 5 * f - 2;
    elseif f <= 0.7,  mu(1) = 1;
    elseif f <= 0.8,  mu(1) = -10 * f + 8;
    else,             mu(1) = 0;
    end
    % S2 - Khai thac
    if f <= 0.2,      mu(2) = 0;
    elseif f <= 0.3,  mu(2) = 10 * f - 2;
    elseif f <= 0.4,  mu(2) = 1;
    elseif f <= 0.6,  mu(2) = -5 * f + 3;
    else,             mu(2) = 0;
    end
    % S3 - Hoi tu
    if f <= 0.1,      mu(3) = 1;
    elseif f <= 0.3,  mu(3) = -5 * f + 1.5;
    else,             mu(3) = 0;
    end
    % S4 - Nhay thoat
    if f <= 0.7,      mu(4) = 0;
    elseif f <= 0.9,  mu(4) = 5 * f - 3.5;
    else,             mu(4) = 1;
    end

    cand = find(mu > 0);
    if isempty(cand)
        [~, state] = max(mu);
    elseif numel(cand) == 1
        state = cand;
    elseif any(cand == prev)
        state = prev;
    else
        % Khoang cach vong tren chu trinh 4 trang thai
        dist = min(mod(cand - prev, 4), mod(prev - cand, 4));
        [~, i] = min(dist);
        state = cand(i);
    end
end
