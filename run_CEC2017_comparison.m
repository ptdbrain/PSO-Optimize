function R = run_CEC2017_comparison(opts)
% RUN_CEC2017_COMPARISON  So sanh PSO, APSO, DE, JADE tren bo CEC2017.
%
% Cach dung:
%   run_CEC2017_comparison              % Cau hinh mac dinh (D = 30, 30 lan chay)
%   run_CEC2017_comparison(struct('runs', 2, 'funcs', [1 3]))  % Chay thu nhanh
%
% Cau hinh mac dinh:
%   D = 30, NP = 30, MaxFEs = 10000*D = 300000, MaxGen = 10000, runs = 30,
%   ham F1, F3-F30 (F2 bi loai khoi CEC2017 chinh thuc), mien [-100, 100]^D.
%
% Ket qua (thu muc results/D30/):
%   raw/Fxx.mat              - ket qua tho tung ham (dung de tiep tuc khi bi ngat)
%   CEC2017_D30_raw.mat      - toan bo ket qua tho
%   CEC2017_D30_compare.xlsx - bang so sanh cac thuat toan
%   figures/*.png            - duong cong hoi tu

    root = fileparts(mfilename('fullpath'));
    cd(root);
    addpath(fullfile(root, 'algorithms'), fullfile(root, 'utils'));

    % ---------------- Cau hinh ----------------
    cfg.D = 30;
    cfg.NP = 30;
    cfg.MaxFEs = 10000 * cfg.D;
    cfg.MaxGen = 10000;
    cfg.runs = 30;
    cfg.lb = -100;
    cfg.ub = 100;
    cfg.funcs = [1, 3:30];
    cfg.algs = {'PSO', 'APSO', 'DE', 'JADE'};
    cfg.ref_alg = 'APSO';                  % Thuat toan tham chieu cho Wilcoxon
    cfg.n_cp = 1000;                       % So diem lay mau duong hoi tu
    cfg.cec_points = [0.01 0.02 0.03 0.05 0.1:0.1:1.0];  % Moc FEs chuan CEC
    cfg.use_parallel = true;               % Dung parfor neu co Parallel Toolbox
    cfg.out_dir = '';

    if nargin > 0
        names = fieldnames(opts);
        for k = 1:numel(names)
            cfg.(names{k}) = opts.(names{k});
        end
        if ~isfield(opts, 'MaxFEs') && isfield(opts, 'D')
            cfg.MaxFEs = 10000 * cfg.D;
        end
    end
    if isempty(cfg.out_dir)
        cfg.out_dir = fullfile(root, 'results', sprintf('D%d', cfg.D));
    end

    if exist('cec17_func', 'file') ~= 3
        error('Khong tim thay cec17_func.mexw64. Hay chay: mex cec17_func.cpp -DWINDOWS');
    end

    raw_dir = fullfile(cfg.out_dir, 'raw');
    if ~exist(raw_dir, 'dir'), mkdir(raw_dir); end

    cp_fes = round(linspace(cfg.MaxFEs / cfg.n_cp, cfg.MaxFEs, cfg.n_cp));
    cec_fes = round(cfg.cec_points * cfg.MaxFEs);

    % ---------------- Song song ----------------
    n_workers = 0;
    if cfg.use_parallel && license('test', 'Distrib_Computing_Toolbox') ...
            && ~isempty(ver('parallel'))
        pool = gcp('nocreate');
        if isempty(pool), pool = parpool; end
        n_workers = pool.NumWorkers;
        % cec17_func doc file theo duong dan tuong doi input_data/
        wait(parfevalOnAll(pool, @cd, 0, root));
        wait(parfevalOnAll(pool, @addpath, 0, fullfile(root, 'algorithms'), fullfile(root, 'utils')));
    end

    fprintf('CEC2017 | D = %d | NP = %d | MaxFEs = %d | MaxGen = %d | runs = %d | workers = %d\n', ...
            cfg.D, cfg.NP, cfg.MaxFEs, cfg.MaxGen, cfg.runs, n_workers);

    algs = cfg.algs;
    nA = numel(algs);
    nF = numel(cfg.funcs);
    nC = numel(cec_fes);

    % ---------------- Chay thi nghiem ----------------
    for i = 1:nF
        fn = cfg.funcs(i);
        raw_file = fullfile(raw_dir, sprintf('F%02d.mat', fn));
        if exist(raw_file, 'file')
            S = load(raw_file, 'cfg_f');
            % Bo qua cac truong khong anh huong ket qua so (thu muc, song song)
            ignore = {'out_dir', 'use_parallel'};
            if isequal(rmfield(S.cfg_f, ignore), rmfield(cfg, ignore))
                fprintf('F%-2d da co ket qua, bo qua.\n', fn);
                continue;
            end
        end

        res = struct();
        for a = 1:nA
            name = algs{a};
            err = zeros(1, cfg.runs);
            X = zeros(cfg.runs, cfg.D);
            curve = zeros(cfg.runs, cfg.n_cp);
            cec_err = zeros(cfg.runs, nC);
            tm = zeros(1, cfg.runs);
            fes = zeros(1, cfg.runs);
            gens = zeros(1, cfg.runs);
            t0 = tic;
            parfor (r = 1:cfg.runs, n_workers)
                % Cung seed cho moi thuat toan tren cung ham/lan chay
                o = single_run(name, fn, cfg, 1000 * fn + r, cp_fes, cec_fes); %#ok<PFBNS>
                err(r) = o.err;
                X(r, :) = o.x;
                curve(r, :) = o.curve;
                cec_err(r, :) = o.cec_err;
                tm(r) = o.time;
                fes(r) = o.fes;
                gens(r) = o.gens;
            end
            res.(name) = struct('err', err, 'x', X, 'curve', curve, ...
                'cec_err', cec_err, 'time', tm, 'fes', fes, 'gens', gens);
            fprintf('F%-2d %-5s mean = %.4e  std = %.4e  (%.1fs)\n', ...
                    fn, name, mean(err), std(err), toc(t0));
        end
        cfg_f = cfg; %#ok<NASGU>
        save(raw_file, 'res', 'cfg_f', 'cp_fes', 'cec_fes');
    end

    % ---------------- Gop ket qua ----------------
    R.cfg = cfg;
    R.cp_fes = cp_fes;
    R.cec_fes = cec_fes;
    for a = 1:nA
        R.(algs{a}) = struct( ...
            'err', zeros(nF, cfg.runs), ...
            'x', zeros(nF, cfg.runs, cfg.D), ...
            'curve', zeros(nF, cfg.runs, cfg.n_cp), ...
            'cec_err', zeros(nF, cfg.runs, nC), ...
            'time', zeros(nF, cfg.runs), ...
            'fes', zeros(nF, cfg.runs), ...
            'gens', zeros(nF, cfg.runs));
    end
    for i = 1:nF
        S = load(fullfile(raw_dir, sprintf('F%02d.mat', cfg.funcs(i))), 'res');
        for a = 1:nA
            s = S.res.(algs{a});
            R.(algs{a}).err(i, :) = s.err;
            R.(algs{a}).x(i, :, :) = s.x;
            R.(algs{a}).curve(i, :, :) = s.curve;
            R.(algs{a}).cec_err(i, :, :) = s.cec_err;
            R.(algs{a}).time(i, :) = s.time;
            R.(algs{a}).fes(i, :) = s.fes;
            R.(algs{a}).gens(i, :) = s.gens;
        end
    end

    tag = sprintf('CEC2017_D%d', cfg.D);
    save(fullfile(cfg.out_dir, [tag '_raw.mat']), '-struct', 'R', '-v7.3');
    fprintf('Da luu ket qua tho: %s\n', fullfile(cfg.out_dir, [tag '_raw.mat']));

    xlsx_file = fullfile(cfg.out_dir, [tag '_compare.xlsx']);
    export_excel(R, xlsx_file);
    fprintf('Da xuat bang so sanh: %s\n', xlsx_file);

    plot_convergence(R, fullfile(cfg.out_dir, 'figures'));
    fprintf('Da luu do thi hoi tu: %s\n', fullfile(cfg.out_dir, 'figures'));
end
