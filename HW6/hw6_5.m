%% CVX_HW6_5

clear all;close all;clc;
%% Setting some parameters

P=[13 ,12,-2; 12,17,6; -2,6,12];
q=[-22; -14.5; 13];
r=1;

alpha= 0.01;
beta= 0.5;
tolerance= 1e-6;
u=[2,10,50,150];

%% Program
figure();
for j=1:length(u)
    x=[0;0;0];
    t=1;
    m=6;
    iter=0;
    i=1;
    duality_gap=[];
    iteration=[];
    while(1)
        % Centering Step
        while(1)
            % Compute Newton step and decrement
            [Xnt ,lambda_square, grad_f]= Newton_step(x,P,q,t);
            % Stopping criterion
            if lambda_square/2 <= tolerance
                break;
            end
            
            % line search 
            l=1; % The original 't' in BTLS.We change it to l
            while(1)
                [f_new, f_old] = BTLS(x,P,q,r,t,l,Xnt);
                if(f_new <= (f_old + alpha*l*grad_f.'*Xnt) )
                    if(-1 < x +l*Xnt & x+l*Xnt<1)
                        break;
                    end
                end
                l = beta * l;

            end
            
            % update
            x = x + l*Xnt;
            iter = iter + 1;
        end
        iteration(i) = iter;
        % Stopping criterion
        duality_gap(i) = m/t;
        if(duality_gap(i) <= tolerance)
            break;
        end

        % increase t 
        t = u(j) * t;
        i = i+1;
    end
    
    fprintf("The optimal point with u = %d\n" , u(j));
    disp(x);
    stairs([0,iteration],[duality_gap(1),duality_gap],'LineWidth',2);
    set(gca,'YScale','log');
    hold on;
end

title('Barrier method');
legend('u=2','u=10','u=50','u=150');
xlabel("Newton iteration");
ylabel("duality gap");
grid on;


function [Xnt, lamda_square, grad_f]=Newton_step(x,P,q,t)
    grad_tf0 = t*P*x + t*q;
    grad2_tf0 = t*P;
    grad_phi = 1./(1-x)-1./(1+x);
    grad2_phi = diag(1./((1-x).^2) + 1./((1+x).^2));
    grad_f = grad_tf0 + grad_phi;
    grad2_f = grad2_tf0 + grad2_phi;
    Xnt = -inv(grad2_f) * grad_f;
    lamda_square = Xnt.' * grad2_f * Xnt;
end

function [f_new, f_old] = BTLS(x,P,q,r,t,l,Xnt)
    f_old = (t/2)*(x.')*P*x + t*q.'*x + t*r -log(1-x(1))-log(1-x(2))-...
    log(1-x(3))-log(1+x(1))-log(1+x(2))-log(1+x(3));
    x_new = x + l*Xnt;
    f_new = (t/2)*(x_new.')*P*x_new + t*q.'*x_new + t*r -...
    log(1-x_new(1))-log(1-x_new(2))-log(1-x_new(3))-log(1+x_new(1))-log(1+x_new(2))-log(1+x_new(3));
end

