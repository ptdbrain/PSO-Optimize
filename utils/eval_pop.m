function [fit, FEs] = eval_pop(fhd, X, func_num, FEs, MaxFEs)
% EVAL_POP  Danh gia quan the (NP x D) va dem so lan danh gia ham.
% Khong bao gio vuot qua MaxFEs: neu ngan sach con lai < NP thi chi danh gia
% cac ca the dau tien, cac ca the con lai nhan gia tri Inf (khong duoc chon).

    NP = size(X, 1);
    k = min(NP, MaxFEs - FEs);
    fit = inf(NP, 1);
    if k > 0
        fit(1:k) = feval(fhd, X(1:k, :)', func_num)';
        FEs = FEs + k;
    end
end
