clear all; close all; clc;
A0 = [10 8 12 15 15; 8 14 8 7 9; 12 8 10 13 9; 15 7 13 4 10; 15 9 9 10 4];
A1 = [12 11 14 10 3; 11 14 10 14 6; 14 10 16 18 4; 10 14 18 18 8; 3 6 4 8 8];
A2 = [4 13 12 16 6; 13 4 14 9 15; 12 14 6 5 5; 16 9 5 2 6; 6 15 5 6 8];

cvx_begin sdp
    variables x1 x2 lamda1 lamda2
    minimize( lamda1-lamda2 )
    subject to
        A = A0 + x1*A1 + x2*A2;
        A <= lamda1*eye(5);
        A >= lamda2*eye(5);
cvx_end

fprintf("x1 = %f, x2 = %f\n",x1,x2);
