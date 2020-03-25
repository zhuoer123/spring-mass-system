clc;clear;
close all;
%%
%%%%%% Spring-mass system
k1=2.5; k2=3.7; %%%N/m
m=1.5;  % kg
h=5.4;  % m
w =sqrt((k1+k2)/m);
A=[0 1;w^2 0];
R=[0.005^2  0;
   0       0.05^2];
Measurement=load('shujv.txt');
time=Measurement(:,1); rou=Measurement(:,2); rou_dot=Measurement(:,3);
X_star=[4.0;0.2];  P=[10000 0;0 100];
deltaX=[0;0];
X=[3.0;0.0];
S1=zeros(2,2); S2=zeros(2,1);S3=0;S4=0;
I=eye(2);
%%
%%%%加权最小二乘法求解
i=0;
 for t=0:length(time)-1
    X(t+1,1)=X(1)*cos(w*t)+X(2)/w*sin(w*t);
    X(t+1,2)=X(2)*cos(w*t)-X(1)*w*sin(w*t);
    
    phi=[cos(w*t)     sin(w*t)/w;
         -w*sin(w*t)  cos(w*t)];   %%% STM
    dX=[X(t+1,1)-X_star(1);X(t+1,2)-X_star(2)];
    dX=phi*dX;
    H_piao=H(X(t+1,1),X(t+1,2),rou(t+1));
    H1=H_piao*phi;
%     y=H_piao*dX;
    y=[rou(t+1,1)-sqrt(X(t+1,1)^2+h^2);rou_dot(t+1,1)-X(t+1)*X(t+2)/rou(t+1,1)];
    S1=S1+H1'*R^-1*H1;
    S2=S2+H1'*R^-1*y; 
 end
 while i<=1000
    deltaX_estimate=inv(S1+P^-1)*(S2+P^-1*deltaX);
    
    X_estimate=X_star+deltaX_estimate;
    deltaX=deltaX-deltaX_estimate;
    if abs(X_estimate-X(1,:)')<=10^-4
        break;
    else 
    end
    deltax(:,i+1)=deltaX_estimate; 
    i=i+1;
 end
 for t=0:length(time)-1
     W=R^-1;
    X(t+1,1)=X(1)*cos(w*t)+X(2)/w*sin(w*t);
    X(t+1,2)=X(2)*cos(w*t)-X(1)*w*sin(w*t);
    
    phi=[cos(w*t)     sin(w*t)/w;
         -w*sin(w*t)  cos(w*t)];   %%% STM

    dX=[X(t+1,1)-X_star(1);X(t+1,2)-X_star(2)];
    dX=phi*dX;
    H_piao=H(X(t+1,1),X(t+1,2),rou(t+1));
    H1=H_piao*phi;
    y=[rou(t+1,1)-sqrt(X(t+1,1)^2+h^2);rou_dot(t+1,1)-X(t+1)*X(t+2)/rou(t+1,1)];
    e=y-H1*X_estimate;
    S3=S3+e'*W*e;
%     S3=S3+e(1)'*W(1,1)*e(1);
%     S4=S4+e(2)'*W(2,2)*e(2);
 end
RMS1=rms(e(1));
RMS2=rms(e(2));
%%
%%%Kalman Filter 求解
P=diag([1000;10].^2);
Q=diag([0.01^2;0.1^2]);
%     选取初值
    Xstar=X_star*ones(1,11);
    deX=X-Xstar';
for ii=0:length(time)-1
% ii=1;
   for t=0:length(time)-1
       phi=[cos(w*t)     sin(w*t)/w;
           -w*sin(w*t)  cos(w*t)];   %%% STM
       dXkf=deX(ii+1,:)';
       dXkf=phi*dXkf;
       H_piao=H(X(t+1,1),X(t+1,2),rou(t+1));
       H1=H_piao*phi;
       %%%量测值
       P=phi*P*phi'+Q;
       z=[rou(t+1,1)-sqrt(X(t+1,1)^2+h^2);rou_dot(t+1,1)-X(t+1)*X(t+2)/rou(t+1,1)];
       K=P*H_piao'*inv(H_piao*P*H_piao'+R);
       dXkf=dXkf+K*(z-H_piao*dXkf);
       P=(I-K*H_piao)*P;
       dX_stor(:,t+1)=dXkf;
       
   end
   flag1=find(min(dX_stor(1,:)));
       DX=dX_stor(:,flag1);
    dX_sum(ii+1,:)=DX;
end
flag=find(min(dX_sum(:,1)));
dX_best=dX_sum(flag,:);
X_kf=X_star+dX_best';





