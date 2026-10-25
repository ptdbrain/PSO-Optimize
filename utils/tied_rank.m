function [r, tiesum] = tied_rank(v)
% TIED_RANK  Xep hang tang dan, cac gia tri bang nhau nhan hang trung binh.
% tiesum = sum(t^3 - t) tren cac nhom bang nhau (dung cho hieu chinh).

    v = v(:);
    n = numel(v);
    [s, ord] = sort(v);
    r = zeros(n, 1);
    tiesum = 0;
    i = 1;
    while i <= n
        j = i;
        while j < n && s(j + 1) == s(i)
            j = j + 1;
        end
        r(ord(i:j)) = (i + j) / 2;
        t = j - i + 1;
        tiesum = tiesum + t ^ 3 - t;
        i = j + 1;
    end
end
