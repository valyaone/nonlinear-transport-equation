clear; clc;
close all;

Uleft = @(t) exp(-t);
Uinit = @(x) (-x + 1);

N = 50;
M = 50;
T = 0.3; %будет косяк если ставить больше 0,17 (либо не будет выводить, либо большие значения слишком)
a = 0;
b = 1;
h = (b - a)/( N - 1);
tau = T/(M-1);

U=zeros(M,N);
tn=0:tau:T;
xn=0:h:(b - a);
   
%начальное условие
for n=1:N 
        U(1,n) = Uinit( xn(n) );
end

%граничное условие
for m=1:M
    U(m,1) = Uleft( tn(m) );
end
    
method = 2; % 1 - ERK1,  2 - CROS1, 3 - CROS2 

if ( method == 1 ) % ERK1
    for k=1:M-1
       w(k,2:N) = F_pr_ch( U(k,2:N), tn(k) + tau/2, h);
       U(k+1,2:N) = U(k,2:N) + tau*w(k,2:N);
    end
end

if (method == 2 ) %CROS1
    A = (1+1i)/2;
    for k=1:M-1
        F = F_pr_ch( U(k,2:N), tn(k) + tau/2, h);
        Fu = yakobian( U(k,2:N), tn(k), h);
        w = F / ( eye(N-1) - A*tau*Fu );
        U(k+1,2:N) = U(k,2:N) + tau*real(w);
    end
end
   
if (method == 3 ) %CROS2
    b1 = 0.1941430241155180-0.2246898944678803i;
    b2 = 0.8058569758844820-0.8870089521907592i;
    c21 = 0.2554708972958462-0.2026195833570109i;
    a21 = 0.5617645150714754-1.148223341045841i;
    a11 = 0.1+sqrt(11)/30*1i;
    a22 = 0.2+0.1i;
    for k=1:M-1
    F1  = F_pr_ch( U(k,2:N), tn(k) + tau/2, h);
    Fu1 = yakobian( U(k,2:N), tn(k), h);
    w1  = F1 / ( eye(N-1) - a11*tau*Fu1 ); 
    
    F2  = F_pr_ch( U(k,2:N) + tau*real( c21*w1 ), tn(k) + tau/2, h);
    Fu2 = yakobian( U(k,2:N) + tau*real( a21*w1 ), tn(k), h);
    w2  = F2 / ( eye(N-1) - a22*tau*Fu2 );
    
    U(k+1,2:N) = U(k,2:N) + tau*real( b1*w1 + b2*w2);
    end
end
    
%Построение 3d графика
fig1 = figure(1);
surf(xn',tn,U);
view([120 45]);
xlabel('X');
ylabel('T');
zlabel('U');
title('Решение уравнения');

fig = figure(2);
% создание первого пустого кадра
set(fig,'Position',[350,200,700,300]);
axes('xlim',[0 1],'ylim',[0 3],'NextPlot','add','Parent',fig);
xlabel('X');
ylabel('U');
title('Решение уравнения');
frame = getframe(fig);
[im,map] = rgb2ind(frame.cdata,4);
imwrite(im,map,'animation2.gif','DelayTime',0,'Loopcount',0);

% цикл анимации
for i=1:M 
    plot(xn',U(i,:),'-');
    %plot(X2(i),Y2(i),'r.');
    
    frame = getframe(fig);
    [im,map] = rgb2ind(frame.cdata,4);
    imwrite(im,map,'animation2.gif','DelayTime',0.1,'WriteMode','Append');
end;

fig3 = figure(3);
%анимация
for i= 1:1:M+1
    plot(xn',U(i,:),'r-' ,'LineWidth' , 4 ) ;
    hold on ;
    %plot ( x(1:1:N+1) ,U(1 : 1 :N+1,i) ,'-', 'LineWidth' , 3 ) ;
    hold on ;
    title('Решение задачи');
    xlabel ('x' ) ;
    ylabel ('U' ) ;
    hold off ;
    drawnow ; 
    pause ( 0.25 ) ;
end