function R = run_CEC2017_comparison(opts)
% RUN_CEC2017_COMPARISON  So sanh cac thuat toan toi uu tren bo CEC2017.
%
% Cach dung:
%   run_CEC2017_comparison                  % Cau hinh mac dinh
%   run_CEC2017_comparison(struct('runs', 3, 'funcs', [1 3], 'out_dir', 'results/test'))
%
% Cau hinh mac dinh (theo chuan CEC2017):
%   D = 30, MaxFEs = 10000*D = 300000, runs = 51, ham F1, F3-F30
%   (F2 bi loai), mien [-100, 100]^D. Chi dung theo so FEs, khong gioi han
%   so the he. Moi thuat toan dung kich thuoc quan the theo tac gia.
%
% Ket qua (thu muc results/D30/):
%   raw/<ALG>/Fxx.mat        - ket qua tho tung thuat toan, tung ham
%   CEC2017_D30_raw.mat      - toan bo ket qua tho
%   CEC2017_D30_compare.xlsx - bang so sanh cac thuat toan
%   figures/*.png            - duong cong hoi tu
%
% Ket qua duoc luu rieng cho tung cap (thuat toan, ham): khi them thuat toan
% moi chi thuat toan do phai chay; neu sua code mot thuat toan, xoa thu muc
% raw/<ALG> cua no de chay lai.

    root = fileparts(mfilename('fullpath'));
    cec_dir = fullfile(root, 'cec2017');
    paths = {fullfile(root, 'algorithms'), fullfile(root, 'algorithms', 'legacy'), ...
             fullfile(root, 'utils'), cec_dir};
    addpath(paths{:});
    % cec17_func doc input_data/ theo duong dan tuong doi
    old_dir = cd(cec_dir);
    restore_dir = onCleanup(@() cd(old_dir));

    % ---------------- Cau hinh ----------------
    cfg.D = 30;
    cfg.MaxFEs = 10000 * cfg.D;
    cfg.runs = 51;
    cfg.lb = -100;
    cfg.ub = 100;
    cfg.funcs = [1, 3:30];
    cfg.algs = {'APSO', 'CLPSO', 'INFO', 'LSHADE'};
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
    elseif ~is_absolute(cfg.out_dir)
        cfg.out_dir = fullfile(root, cfg.out_dir);
    end

    if exist('cec17_func', 'file') ~= 3
        error('Khong tim thay cec17_func.mexw64. Trong thu muc cec2017 hay chay: mex cec17_func.cpp -DWINDOWS');
    end

    raw_dir = fullfile(cfg.out_dir, 'raw');
    cp_fes = round(linspace(cfg.MaxFEs / cfg.n_cp, cfg.MaxFEs, cfg.n_cp));
    cec_fes = round(cfg.cec_points * cfg.MaxFEs);

    % Cac truong cau hinh quyet dinh ket qua so (dung de kiem tra file cu)
    key = struct('D', cfg.D, 'MaxFEs', cfg.MaxFEs, 'runs', cfg.runs, ...
                 'lb', cfg.lb, 'ub', cfg.ub, 'n_cp', cfg.n_cp, 'cec_points', cfg.cec_points);

    % ---------------- Song song ----------------
    n_workers = 0;
    if cfg.use_parallel && license('test', 'Distrib_Computing_Toolbox') ...
            && ~isempty(ver('parallel'))
        pool = gcp('nocreate');
        if isempty(pool), pool = parpool; end
        n_workers = pool.NumWorkers;
        wait(parfevalOnAll(pool, @addpath, 0, paths{:}));
        wait(parfevalOnAll(pool, @cd, 0, cec_dir));
    end

    fprintf('CEC2017 | D = %d | MaxFEs = %d | runs = %d | %s | workers = %d\n', ...
            cfg.D, cfg.MaxFEs, cfg.runs, strjoin(cfg.algs, ', '), n_workers);

    algs = cfg.algs;
    nA = numel(algs);
    nF = numel(cfg.funcs);
    nC = numel(cec_fes);

    % ---------------- Chay thi nghiem ----------------
    for i = 1:nF
        fn = cfg.funcs(i);
        for a = 1:nA
            name = algs{a};
            raw_file = fullfile(raw_dir, name, sprintf('F%02d.mat', fn));
            if exist(raw_file, 'file')
                S = load(raw_file, 'key');
                if isequal(S.key, key)
                    fprintf('F%-2d %-6s da co ket qua, bo qua.\n', fn, name);
                    continue;
                end
            end

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
            res = struct('err', err, 'x', X, 'curve', curve, 'cec_err', cec_err, ...
                         'time', tm, 'fes', fes, 'gens', gens); %#ok<NASGU>
            if ~exist(fileparts(raw_file), 'dir'), mkdir(fileparts(raw_file)); end
            save(raw_file, 'res', 'key');
            fprintf('F%-2d %-6s mean = %.4e  std = %.4e  (%.1fs)\n', ...
                    fn, name, mean(err), std(err), toc(t0));
        end
    end

    % ---------------- Gop ket qua ----------------
    R.cfg = cfg;
    R.cp_fes = cp_fes;
    R.cec_fes = cec_fes;
    for a = 1:nA
        A = struct( ...
            'err', zeros(nF, cfg.runs), ...
            'x', zeros(nF, cfg.runs, cfg.D), ...
            'curve', zeros(nF, cfg.runs, cfg.n_cp), ...
            'cec_err', zeros(nF, cfg.runs, nC), ...
            'time', zeros(nF, cfg.runs), ...
            'fes', zeros(nF, cfg.runs), ...
            'gens', zeros(nF, cfg.runs));
        for i = 1:nF
            S = load(fullfile(raw_dir, algs{a}, sprintf('F%02d.mat', cfg.funcs(i))), 'res');
            s = S.res;
            A.err(i, :) = s.err;
            A.x(i, :, :) = s.x;
            A.curve(i, :, :) = s.curve;
            A.cec_err(i, :, :) = s.cec_err;
            A.time(i, :) = s.time;
            A.fes(i, :) = s.fes;
            A.gens(i, :) = s.gens;
        end
        R.(algs{a}) = A;
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

function tf = is_absolute(p)
    tf = numel(p) >= 2 && (p(2) == ':' || p(1) == '\' || p(1) == '/');
end
