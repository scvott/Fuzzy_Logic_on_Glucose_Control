clc
clear all
close all

%parameter of membership function for glucose
	gb = 70;a=-25; b=10; step=0.1; t=a:step:b;
	for n = 1:length(t)
	A = [1,(-15-t(n))/(-15-(-20))];
    B = [(t(n)-(-20))/((-15)-(-20)),(-10-t(n))/(-10-(-15))];
	C = [(t(n)-(-15))/((-10)-(-15)),(-5-t(n))/(-5-(-10))];
	D = [(t(n)-(-10))/(-7-(-10)),(-2-t(n))/(-2-( -7))];
	E = [(t(n)-(-5))/(-2-(-5)),(0-t(n))/(0+2)];
	F = [(t(n)-(-2))/(0-(-2)),(2-t(n))/(2-0)];
	G = [(t(n)-0)/(5-0),1];
	
	f1(n) = max(min(A),0);f2(n) = max(min(B),0);f3(n) = max(min(C),0);f4(n) = max(min(D),0);f5(n) = max(min(E),0);f6(n) = max(min(F),0);f7(n) = max(min(G),0);
	end
	figure(1)
    %title("Membership Function");
    xlabel('glucose error'); 
    ylabel('degree of membership');
	hold on
	plot(t,f1(1:length(t)),'LineWidth',1.5);
	plot(t,f2(1:length(t)),'LineWidth',1.5);
	plot(t,f3(1:length(t)),'LineWidth',1.5);
	plot(t,f4(1:length(t)),'LineWidth',1.5);
	plot(t,f5(1:length(t)),'LineWidth',1.5);
	plot(t,f6(1:length(t)),'LineWidth',1.5);
	plot(t,f7(1:length(t)),'LineWidth',1.5);
	hold off
    legend('NVL','NL','NB','NM','NS','NOM','P');  