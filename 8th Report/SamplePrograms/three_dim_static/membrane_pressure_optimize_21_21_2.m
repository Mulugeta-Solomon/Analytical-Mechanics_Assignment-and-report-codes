% 平衡式と制約式をまとめた等式制約とする
% 外力 fn は圧力による力であり，定ベクトルではないので
% 等式制約を nonlcon に記述する
% g, cm, sec

dt0 = datetime;
addpath('../three_dim_fea');