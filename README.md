# PSO-Optimize

So sánh **APSO**, **CLPSO**, **INFO** và **L-SHADE** trên bộ hàm benchmark **CEC2017** (MATLAB).

## Cấu trúc thư mục

```
├── run_CEC2017_comparison.m   Script chính: chạy thí nghiệm, lưu .mat, xuất .xlsx, vẽ hình
├── algorithms/                Các thuật toán (chỉ tính toán, không ghi file)
│   ├── alg_APSO.m             Adaptive PSO (Zhan et al., 2009)
│   ├── alg_CLPSO.m            Comprehensive Learning PSO (Liang et al., 2006)
│   ├── alg_INFO.m             INFO - weighted mean of vectors (Ahmadianfar et al., 2022)
│   ├── alg_LSHADE.m           L-SHADE (Tanabe & Fukunaga, 2014)
│   └── legacy/                Thuật toán cũ, không chạy mặc định
│       ├── alg_PSO.m          PSO trọng số quán tính giảm tuyến tính
│       ├── alg_DE.m           DE/rand/1/bin
│       └── alg_JADE.m         JADE
├── utils/                     Hàm hỗ trợ
│   ├── eval_pop.m             Đánh giá quần thể, đếm FEs (không vượt MaxFEs)
│   ├── repair_bounds.m        Xử lý nghiệm vượt biên (L-SHADE, DE, JADE)
│   ├── single_run.m           Chạy một lần một thuật toán trên một hàm
│   ├── export_excel.m         Xuất bảng so sánh ra .xlsx
│   ├── plot_convergence.m     Vẽ đường cong hội tụ
│   ├── rank_sum_test.m        Kiểm định Wilcoxon rank-sum
│   └── tied_rank.m            Xếp hạng có xử lý giá trị bằng nhau
├── cec2017/                   Bộ hàm CEC2017
│   ├── cec17_func.cpp         Mã nguồn C++
│   ├── cec17_func.mexw64      Bản biên dịch cho MATLAB trên Windows
│   └── input_data/            Dữ liệu dịch chuyển, xoay, hoán vị
└── results/D30/               Kết quả: bảng .xlsx và hình hội tụ
```

## Cấu hình thí nghiệm

Theo chuẩn CEC2017:

| Tham số | Giá trị |
|---|---|
| Số chiều D | 30 |
| MaxFEs | 10000 × D = 300000 (điều kiện dừng duy nhất) |
| Số lần chạy | 51 |
| Hàm | F1, F3–F30 (F2 bị loại theo chuẩn CEC2017) |
| Miền tìm kiếm | [-100, 100]^D |
| Sai số | f(x) − f(x*), giá trị < 1e-8 coi là 0 |

Mỗi thuật toán dùng kích thước quần thể và tham số theo tác giả:

| Thuật toán | Kích thước quần thể | Tham số chính | Nguồn cài đặt |
|---|---|---|---|
| APSO | 20 | w(f) ∈ [0.4, 0.9], c1, c2 ∈ [1.5, 2.5], δ ∈ [0.05, 0.1], σ: 1.0 → 0.1 | Tự cài theo bài báo (không có code công khai) |
| CLPSO | 40 | w: 0.9 → 0.2, c = 1.49445, Pc ∈ [0, 0.5], chọn lại hạt mẫu sau 5 thế hệ | Chuyển từ code MATLAB gốc `CLPSO_new_func.m` |
| INFO | 30 | α = 2·exp(−4·FEs/MaxFEs) | Chuyển từ code MATLAB gốc `INFO.m` |
| L-SHADE | 18·D → 4 (giảm tuyến tính) | H = 6, p = 0.11, \|A\| = 2.6·NP | Chuyển từ code C++ gốc `lshade.cc` (CEC2014) |

Với các thuật toán chuyển từ code gốc, chỉ phần giao tiếp được thay đổi (đếm FEs, ghi đường hội tụ, dừng theo MaxFEs); phần thuật toán giữ nguyên, kể cả những chỗ code khác bài báo. Nguồn code:

- CLPSO: kho [P-N-Suganthan/CODES](https://github.com/P-N-Suganthan/CODES), file `2006-IEEE-TEC-CLPSO.zip`.
- INFO: [trang dự án INFO](https://mdm.wzu.edu.cn/pages/algorithms/INFO/INFO.html), file `Matlab-codes-of-INFO-optimization-algorithm.zip`.
- L-SHADE: kho [P-N-Suganthan/CEC2014](https://github.com/P-N-Suganthan/CEC2014), file `Top-Methods-Part-A.rar` → `LSHADE_CEC2014.zip`.

## Cách chạy

Yêu cầu: MATLAB (đã thử trên R2019b). Parallel Computing Toolbox là tùy chọn; nếu có, các lần chạy sẽ chạy song song.

```matlab
run_CEC2017_comparison              % chạy toàn bộ thí nghiệm
```

Chạy thử nhanh với cấu hình nhỏ:

```matlab
run_CEC2017_comparison(struct('funcs', [1 3], 'runs', 3, 'out_dir', 'results/test'))
```

Thêm thuật toán cũ vào so sánh:

```matlab
run_CEC2017_comparison(struct('algs', {{'APSO', 'CLPSO', 'INFO', 'LSHADE', 'JADE'}}))
```

Chạy riêng một thuật toán (kết quả chỉ trả về Workspace, không ghi file):

```matlab
addpath('algorithms', 'algorithms/legacy', 'utils', 'cec2017')
cd cec2017                          % bắt buộc: cec17_func đọc input_data/ theo đường dẫn tương đối
fn = 5;
[best, x, curve, fes] = alg_LSHADE(@cec17_func, fn, 30, 300000, -100, 100);
err = best - 100 * fn
```

Nếu `cec17_func.mexw64` không chạy được (hệ điều hành khác Windows hoặc khác phiên bản MATLAB), biên dịch lại trong thư mục `cec2017`:

```matlab
mex cec17_func.cpp -DWINDOWS
```

## Kết quả

Được ghi vào `results/D30/`:

- `CEC2017_D30_compare.xlsx`: thống kê Mean/Std/Best/Worst/Median, kiểm định Wilcoxon (so với APSO), xếp hạng Friedman, sai số tại các mốc FEs chuẩn CEC, thời gian chạy.
- `figures/`: đường cong hội tụ của từng hàm và hình tổng hợp `All_functions.png`.
- `CEC2017_D30_raw.mat`, `raw/<thuật toán>/Fxx.mat`: kết quả thô. Các file này **không được đưa lên git**, chạy lại script để sinh ra.

Kết quả thô được lưu riêng cho từng cặp (thuật toán, hàm). Khi chạy lại, các cặp đã có kết quả với cùng cấu hình sẽ được bỏ qua. Vì vậy:

- Thêm thuật toán mới: chỉ thuật toán đó phải chạy.
- Thí nghiệm bị dừng giữa chừng: chạy lại sẽ tiếp tục từ chỗ dừng.
- Sửa code một thuật toán: xóa thư mục `results/D30/raw/<thuật toán>` để chạy lại thuật toán đó.

## Tài liệu tham khảo

- Z.-H. Zhan, J. Zhang, Y. Li, H. S.-H. Chung, "Adaptive Particle Swarm Optimization," *IEEE Trans. SMC-B*, 39(6), 1362–1381, 2009.
- J. J. Liang, A. K. Qin, P. N. Suganthan, S. Baskar, "Comprehensive Learning Particle Swarm Optimizer for Global Optimization of Multimodal Functions," *IEEE Trans. Evol. Comput.*, 10(3), 281–295, 2006.
- I. Ahmadianfar, A. A. Heidari, S. Noshadian, H. Chen, A. H. Gandomi, "INFO: An efficient optimization algorithm based on weighted mean of vectors," *Expert Systems with Applications*, 195, 116516, 2022.
- R. Tanabe, A. S. Fukunaga, "Improving the Search Performance of SHADE Using Linear Population Size Reduction," *Proc. IEEE CEC*, 1658–1665, 2014.
- Y. Shi, R. C. Eberhart, "Empirical study of particle swarm optimization," *Proc. CEC*, 1999.
- R. Storn, K. Price, "Differential Evolution – A Simple and Efficient Heuristic for Global Optimization over Continuous Spaces," *J. Global Optimization*, 11, 341–359, 1997.
- J. Zhang, A. C. Sanderson, "JADE: Adaptive Differential Evolution With Optional External Archive," *IEEE Trans. Evol. Comput.*, 13(5), 945–958, 2009.
- N. H. Awad, M. Z. Ali, J. J. Liang, B. Y. Qu, P. N. Suganthan, "Problem Definitions and Evaluation Criteria for the CEC 2017 Special Session and Competition on Single Objective Real-Parameter Numerical Optimization," Technical Report, 2016. Mã nguồn `cec17_func.cpp` và `input_data/` thuộc về các tác giả này.
