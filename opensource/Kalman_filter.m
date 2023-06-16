function [x_est,p_est,y] = Kalman_filter(st,sf,sl,Ts,x_est,p_est,z)

A = [1,Ts,0;
    0,1,0;
    0,0,1;];

R = diag([st,sf,sl]);

Q = [st*Ts+1/3*sf*power(Ts,3),1/2*sf*power(Ts,2),0;
    1/2*sf*power(Ts,2),sf*Ts,0;
    0,0,sl*Ts];

H = eye(3);

x_prd = A * x_est;
p_prd = A * p_est * A' + Q;
K = (p_prd*H')/(H * p_prd * H' + R);
x_est = x_prd + K * (z - H * x_prd);
p_est = p_prd - K * H * p_prd;
y = H * x_est;
end