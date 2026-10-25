function out = single_run(alg_name, func_num, cfg, seed, cp_fes, cec_fes)
% SINGLE_RUN  Chay mot lan doc lap cua mot thuat toan tren mot ham CEC2017.
% Tra ve sai so cuoi, nghiem tot nhat, duong hoi tu (lay mau tai cac moc FEs)
% va thoi gian chay.

    rng(seed, 'twister');
    alg = str2func(['alg_' alg_name]);
    fstar = 100 * func_num;

    t0 = tic;
    [val, x, cv, cf] = alg(@cec17_func, func_num, cfg.D, cfg.NP, ...
                           cfg.MaxFEs, cfg.MaxGen, cfg.lb, cfg.ub);
    out.time = toc(t0);

    out.err = clip_err(val - fstar);
    out.x = x(:)';
    out.fes = cf(end);
    out.gens = numel(cf);
    out.curve = clip_err(sample_curve(cv - fstar, cf, cp_fes));
    out.cec_err = clip_err(sample_curve(cv - fstar, cf, cec_fes));
end

function e = clip_err(e)
% Chuan CEC: sai so < 1e-8 duoc coi la 0
    e(e < 1e-8) = 0;
end

function y = sample_curve(cv, cf, points)
% Gia tri tot nhat tai moi moc FEs (bac thang: gia tri cua the he gan nhat
% co so FEs <= moc).
    y = zeros(1, numel(points));
    for k = 1:numel(points)
        idx = find(cf <= points(k), 1, 'last');
        if isempty(idx), idx = 1; end
        y(k) = cv(idx);
    end
end
