function V = repair_bounds(V, X, LB, UB)
% REPAIR_BOUNDS  Dua vector dot bien vi pham bien ve trung diem giua ca the
% cha va bien (cach xu ly cua JADE, Zhang & Sanderson 2009).
    low = V < LB;
    V(low) = (LB(low) + X(low)) / 2;
    high = V > UB;
    V(high) = (UB(high) + X(high)) / 2;
end
