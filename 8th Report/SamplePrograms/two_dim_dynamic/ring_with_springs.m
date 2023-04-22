% g, cm, sec

addpath('../two_dim_fea');

period01 = 0.2; period12 = 1.8; period23 = 0.1; period34 = 4.8; period45 = 0.1; period56 = 5.0;
% time1=0.2; time2 = 2.0; time3 = 2.1; time4 = 6.6; time5 = 6.7; time6 = 12.0;
time1 =  0 + period01; time2 = time1 + period12; time3 = time2 + period23; time4 = time3 + period34; time5 = time4 + period45; time6 = time5 + period56;
tinterval = 0.1;
vinterval = 0.01;
[ 0, time1, time2, time3, time4, time5, time6 ]