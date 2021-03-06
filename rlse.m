function rlse()
clc
clear
xo=[2 1]';%The true value of unknown coefficient vector

temperature_celsius = [13, 14, 17, 18, 19, 15, 13, 31, 32, 29, 27];
temperature_farenheit = [55, 58, 63, 65, 66, 59, 56, 87, 90, 85, 81];

fprintf('\n\n\temperature_celsius = %f\n\t\t\t %f\n', length(temperature_celsius));
fprintf('\n\n\temperature_farenheit = %f\n\t\t\t %f\n', length(temperature_farenheit));

% (0.01*k)*x1+x2=b

NA=2; % length(temperature_celsius);
x=zeros(NA,1);

for i=1:3
    P=(10^i)*eye(NA,NA);

    fprintf('\n\n\NA = %f\n\t\t\t %f\n', length(NA));

    %Jei k iki 10, o P=1*..., tai gaunasi blogas tikslumas, jei P=10000*...
    %tikslumas gereja...
    xx=[]
    for k=1:length(temperature_celsius), %padidinus iteraciju skaiciu, pvz. iki 1000 sprendinys arteja prie xo=[2 1]
        A(k,:)=[temperature_farenheit(k) 1];
        b(k,:)=temperature_celsius(k); %A(k,:)*xo+0.2*randn;
        [x,K,P]=rlse_online(A(k,:),b(k,:),x,P);
        xx=[xx x];
    end

    fprintf('\n\n\nx ivertis taikant RLSE = %f\n\t\t\t %f\n',x(1),x(2));

    x1=pinv(A)*b; %x1=A\b;
    fprintf('Palyginimui pinv(A)*b = %f\n\t\t  %f\n',x1);

	% jeigu norime atvaizduoti viename grafike ->
%     plot(xx');
%     grid;
%     hold on;

    % atvaizduojame parametru iverciu priklausomybe i grafika pagal
    % stiprinima
    subplot(2,1,1);
    plot(xx(1,:));
    title('Koeficientas k1');
    legend('stiprinimas = 10', 'stiprinimas = 100', 'stiprinimas = 1000');
    grid;
    hold on;
    subplot(2,1,2);
    plot(xx(2,:));
    title('Koeficientas k2');
    legend('stiprinimas = 10', 'stiprinimas = 100', 'stiprinimas = 1000');
    grid;
    hold on;

end

hold off;


function[x,K,P]=rlse_online(aT_k1,b_k1,x,P)
K=P*aT_k1'/(aT_k1*P*aT_k1'+1);%Eq.(18) %nuimti ; pavaizduoti K matricos verte, kuri mazejai, kai iteraciju daugeja
x=x+K*(b_k1-aT_k1*x);%Eq.(17)
P=P-K*aT_k1*P;%Eq.(19)