function p = rank_sum_test(x, y)
% RANK_SUM_TEST  Kiem dinh Wilcoxon rank-sum (Mann-Whitney U) hai phia,
% xap xi chuan co hieu chinh hang bang nhau va hieu chinh lien tuc.
% Khong can Statistics Toolbox.

    x = x(:); y = y(:);
    nx = numel(x); ny = numel(y); N = nx + ny;
    [r, tiesum] = tied_rank([x; y]);

    W = sum(r(1:nx));
    mu = nx * (N + 1) / 2;
    sigma2 = nx * ny / 12 * ((N + 1) - tiesum / (N * (N - 1)));
    if sigma2 <= 0
        p = 1;      % Tat ca gia tri bang nhau
        return;
    end
    z = (W - mu - 0.5 * sign(W - mu)) / sqrt(sigma2);
    p = min(1, erfc(abs(z) / sqrt(2)));
end
