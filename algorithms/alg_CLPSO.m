function [gbestval, gbest, curve_val, curve_fes] = alg_CLPSO(fhd, func_num, D, MaxFEs, lb, ub)
% ALG_CLPSO  PSO hoc toan dien (Comprehensive Learning PSO)
% J. J. Liang, A. K. Qin, P. N. Suganthan, S. Baskar, "Comprehensive Learning
% Particle Swarm Optimizer for Global Optimization of Multimodal Functions",
% IEEE Trans. Evol. Comput., 10(3), 281-295, 2006.
%
% Chuyen the tu code goc CLPSO_new_func.m (kho CODES cua P. N. Suganthan,
% 2006-IEEE-TEC-CLPSO.zip). Chi thay doi phan giao tiep: dem FEs, ghi
% duong hoi tu, dung theo MaxFEs. Thuat toan giu nguyen nhu code goc:
%   - Moi chieu hoc tu pbest cua mot hat mau (dau loai 2 hat) voi xac suat Pc_k
%   - Pc_k = 0.5*(exp(t_k) - 1)/(exp(5) - 1), t_k tu 0 den 5
%   - w giam tuyen tinh 0.9 -> 0.2 theo the he, c = 1.49445
%   - Chon lai hat mau khi pbest khong cai thien 5 the he lien tiep
%   - Cap nhat bat dong bo tung hat; hat ngoai mien khong duoc danh gia
%   NP = 40 cho D = 30 theo bai bao.

    ps = 40;
    me = ceil(MaxFEs / ps);         % Max_Gen = Max_FES / ps nhu TEC_test.m
    cc = 1.49445;
    stay_max = 5;

    t = 5 * (0:ps - 1) / (ps - 1);
    Pc = 0.5 * (exp(t) - exp(t(1))) / (exp(t(ps)) - exp(t(1)));

    lb = repmat(lb, 1, D / numel(lb));
    ub = repmat(ub, 1, D / numel(ub));
    mv = 0.2 * (ub - lb);

    pos = repmat(lb, ps, 1) + repmat(ub - lb, ps, 1) .* rand(ps, D);
    [e, FEs] = eval_pop(fhd, pos, func_num, 0, MaxFEs);
    vel = -repmat(mv, ps, 1) + 2 * repmat(mv, ps, 1) .* rand(ps, D);

    pbest = pos;
    pbestval = e;
    [gbestval, gid] = min(pbestval);
    gbest = pbest(gid, :);

    stay_num = zeros(ps, 1);
    f_pbest = zeros(ps, D);
    for k = 1:ps
        f_pbest(k, :) = assign_exemplar(k, Pc(k), pbestval, D);
    end

    curve_val = zeros(1, MaxFEs);
    curve_fes = zeros(1, MaxFEs);
    gen = 1;
    curve_val(gen) = gbestval;
    curve_fes(gen) = FEs;

    % gen < MaxFEs: chan an toan, tranh lap vo han neu moi hat deu ra ngoai mien
    while FEs < MaxFEs && gen < MaxFEs
        gen = gen + 1;
        w = 0.9 - min(gen, me) * (0.7 / me);

        for k = 1:ps
            if stay_num(k) >= stay_max
                stay_num(k) = 0;
                f_pbest(k, :) = assign_exemplar(k, Pc(k), pbestval, D);
            end

            pbest_f = pbest(sub2ind([ps, D], f_pbest(k, :), 1:D));
            vel(k, :) = w * vel(k, :) + cc * rand(1, D) .* (pbest_f - pos(k, :));
            vel(k, :) = min(max(vel(k, :), -mv), mv);
            pos(k, :) = pos(k, :) + vel(k, :);

            % Chi danh gia hat nam trong mien tim kiem
            if all(pos(k, :) <= ub) && all(pos(k, :) >= lb) && FEs < MaxFEs
                ek = feval(fhd, pos(k, :)', func_num);
                FEs = FEs + 1;
                if pbestval(k) <= ek
                    stay_num(k) = stay_num(k) + 1;
                else
                    pbest(k, :) = pos(k, :);
                    pbestval(k) = ek;
                    if ek < gbestval
                        gbest = pos(k, :);
                        gbestval = ek;
                    end
                end
            end
        end

        curve_val(gen) = gbestval;
        curve_fes(gen) = FEs;
    end

    curve_val = curve_val(1:gen);
    curve_fes = curve_fes(1:gen);
end

function f = assign_exemplar(k, pc, pbestval, D)
% Chon hat mau cho tung chieu cua hat k (nhu code goc)
    ps = numel(pbestval);
    fi1 = ceil(ps * rand(1, D));
    fi2 = ceil(ps * rand(1, D));
    better1 = pbestval(fi1(:))' < pbestval(fi2(:))';
    fi = better1 .* fi1 + ~better1 .* fi2;
    bi = ceil(rand(1, D) - 1 + pc);         % = 1 voi xac suat pc
    if ~any(bi)
        bi(randi(D)) = 1;
    end
    f = bi .* fi + (1 - bi) .* k;
end
