%% CVX_HW6_4(a)

clear all;close all;clc;

%% Setting some parameters

epsilon=1e-8;
start_x=[1;2;3;4;5];
x=start_x;
A=[4,3,2,1,0;1,1,1,1,1];
b=[20;15];
% backtracking line search parameters
alpha=0.01;
beta=0.5;

%% Algorithm

while(1)
    % Compute Newton Step and decrement
    grad_1=log(x)+1;
    grad_2=diag(1./x);
    KKT_sol=inv([grad_2,A.';A,zeros(size(A,1),size(A,1))]) * [(-1)*grad_1;zeros(size(A,1),1)];
    NtStep=KKT_sol(1:length(grad_2(1,:)));
    
    % Stopping criterion
    lamda_square=(NtStep.')*grad_2*(NtStep);
    if(lamda_square/2)<=epsilon
        break;
    end

    % Line search
    StepSize = BTLS(x,grad_1,NtStep,alpha,beta);

    % Update
    x = x + StepSize * NtStep;
end

fprintf("The initial point is : [%s]\n",join(string(start_x)));
fprintf("The optimal point is : [%s]\nThe optimal value is : %.4d\n",join(string(x),','),x.'*log(x));


%% Backtracking line search function

function t= BTLS(x,grad_1,NtStep,alpha,beta)
    t=1;
    f_new=(x+t*NtStep).' * log(x+t*NtStep);
    f_curr=x.' * log(x);
    while f_new > f_curr + alpha * t * (grad_1.' *NtStep);
        t = beta * t;
        f_new = (x+t*NtStep).' * log(x+t*NtStep);
    end
    while sum((x+t*NtStep)>0)~=5
        t= beta * t;
    end
end











