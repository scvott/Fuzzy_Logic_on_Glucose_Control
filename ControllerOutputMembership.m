clear;
a=-2; b=6; step=0.1; t=a:step:b;
for i = 1:length(t)

f1 = [(t(i)-(-2))/((-1)-(-2)),(0-t(i))/(0-(-1)),1];
f2 = [(t(i)-(-1))/(0-(-1)),(0.001-t(i))/(0.001-0),1];
f3 = [(t(i)-0)/(0.001-0),(0.009-t(i))/(0.009-0.001),1];
f4 = [(t(i)-0.001)/(1-0.001),(2-t(i))/(2-1),1];
f5 = [(t(i)-0.009)/(2-0.009),(3-t(i))/(3-2),1];
f6 = [(t(i)-2)/(3-2),(4.5-t(i))/(4.5-3),1];
f7 = [(t(i)-3)/(4.5-3),(6-t(i))/(6-4.5),1];

NS(i) = max(min(f1),0);
Z(i) = max(min(f2),0);
PS(i) = max(min(f3),0);
PM(i) = max(min(f4),0);
PB(i) = max(min(f5),0);
PL(i) = max(min(f6),0);
PVL(i) = max(min(f7),0);
end
figure(1)
    %title("Membership Function");
    xlabel('glucose error'); 
    ylabel('degree of membership');
	hold on
	plot(t,NS(1:length(t)));
	plot(t,Z(1:length(t)));
	plot(t,PS(1:length(t)));
	plot(t,PM(1:length(t)));
	plot(t,PB(1:length(t)));
	plot(t,PL(1:length(t)));
	plot(t,PVL(1:length(t)));
	hold off
    legend('NS','Z','PS','PM','PB','PL','PVL');  

