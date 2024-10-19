%% CVX_HW6_4(b)

clear all;close all;clc;
%% Setting some parameters

epsilon=1e-8;
start_x=[1,5;2,2;3,3;4,4;5,5];
x=start_x;
A=[4,3,2,1,0;1,1,1,1,1];
b=[20;15];
v=zeros(2,1);
% backtracking line search parameters
alpha=0.01;
beta=0.5;
count=0;

%% Algorithm

for i=1:2
    while(1)
        % Compute Newton Step and decrement
        grad_1=log(x(:,i))+1;
        grad_2=diag(1./x(:,i));
        KKT_sol=inv([grad_2,A.';A,zeros(size(A,1),size(A,1))]) * (-1)*[(grad_1 + A.'*v);A*x(:,i) - b];
        size_primal = length(grad_2(1,:));
        size_dual = length(A(:,1));
        NtStep_primal = KKT_sol(1:size_primal);
        NtStep_dual = KKT_sol(size_primal+1:size_primal+size_dual);

        % Line search
        StepSize = BTLS(x,grad_1,A,b,v,i,NtStep_primal,NtStep_dual,alpha,beta);

        % Update
        x(:,i) = x(:,i) + StepSize * NtStep_primal;
        v = v + StepSize * NtStep_dual;

        % Stopping criterion
        r = [(grad_1 + A.'*v); (A * x(:,i) -b)];
        if (A*x(:,i) == b) & (norm(r) <= epsilon)
            break;
        end
    end

end
for i = 1:2
    fprintf("(%d) The initial point is : [%s]\n",i,join(string(start_x(:,i))));
    fprintf("    The optimal point is : [%s]\n    The optimal value is : %.4d\n",join(string(x(:,i)),','),x(:,i).'*log(x(:,i)));
end

%% Backtracking line search function

function t = BTLS(x,grad,A,b,v,i,primal_NewtonStep,dual_NewtonStep,alpha,beta)
    t = 1;
    r = [(grad + A.'*v); A*x(:,i) - b];
    primal_new = (log(x(:,i)+t*primal_NewtonStep) + 1) + A.'*(v+t*dual_NewtonStep);
    while norm([primal_new; A*(x(:,i)+t*primal_NewtonStep)-b]) > (1-alpha*t)*norm(r)
        t = beta * t;
        primal_new = log(x(:,i)+t*primal_NewtonStep) + 1 + A.'*(v+t*dual_NewtonStep);
    end
    %Check
    while sum((x(:,i) + t*primal_NewtonStep) > 0) ~= 5
        t = beta * t;
    end
end








