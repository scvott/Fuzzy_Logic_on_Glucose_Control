clc
clear all
close all
a=0; b=1200; step=0.1; t=a:step:b; %for length
%parameter for people
p1=0; p2=0.0123; p3=5.3*10^(-6);
r=0.005;
n=0.3;
h=78;
gb=70; ib=7;
g(1)=70;
i(1)=7;
x(1)=0;
Da = 0.5; Db = 0.05;%for D(t)
%vlaue for insulin infusion
NS(1) = -0.5;Z(1) = 0;PS(1) = 0.1;PM(1) = 0.4;PB(1) = 1;PL(1) = 2;
pace = 0.01;
%parameter for RFNNI
ata = 0.2;
sigma(1) = 50;    	 %sigma
%weight
w1(1) = 70;w2(1) = 70;w3(1) = 70;w4(1) = 70;w5(1) = 70;w6(1) = 70;w7(1) = 70;w8(1) = 70;w9(1) = 70;
%mid point
m1(1)=50;m2(1)=150;m3(1)=250;m4(1)=-1;m5(1)=0;m6(1)=1;
for j = 1:length(t)
	%eat
	if t(j) >= 800
        d(j)=Da*exp(-Db*(t(j)-800));
    elseif t(j)>= 400
        d(j)=Da*exp(-Db*(t(j)-400));
    else
         d(j)=Da*exp(-Db*t(j));
    end
		%parameter for membership
	xc1 = 70-g(j);				%cintrol input
	xc2 = -(-x(j)*g(j)+d(j));
	%parameter of membership function for glucose
	[line] = LineFunction(0, -20, -15, 1, xc1);
    A = [1,line];
	B = [(xc1-(-15))/((-10)-(-15)),(-5-xc1)/(-5-(-10))];
	C = [(xc1-(-10))/(-7-(-10)),(-2-xc1)/(-2-( -7))];
	D = [(xc1-(-5))/(-2-(-5)),(0-xc1)/(0+2)];
	E = [(xc1-(-2))/(0-(-2)),(2-xc1)/(2-0)];
	F = [(xc1-0)/(5-0),1];
	%membership function for glucose
	function1(j) = max(min(A),0);function2(j) = max(min(B),0);function3(j) = max(min(C),0);function4(j) = max(min(D),0);function5(j) = max(min(E),0);function6(j) = max(min(F),0);
	%parameter of membership function for glucose deviation
	N = [1,(-xc2/20)]; 
    ZERO = [(xc2+1),(1-xc2)];
    P = [xc2/20,1];
	%membership function for glucose deviation
	neg = max(min(N),0);
	zero = max(min(ZERO),0);
	pos = max(min(P),0);
	%rule
	r1=min(function1(j),neg); r2=min(function1(j),zero); r3=min(function1(j),pos);
	r4=min(function2(j),neg); r5=min(function2(j),zero); r6=min(function2(j),pos);
	r7=min(function3(j),neg); r8=min(function3(j),zero); r9=min(function3(j),pos);
	r10=min(function4(j),neg); r11=min(function4(j),zero); r12=min(function4(j),pos);
	r13=min(function5(j),neg); r14=min(function5(j),zero); r15=min(function5(j),pos);
	r16=min(function6(j),neg); r17=min(function6(j),zero); r18=min(function6(j),pos);
	rsum = r1+r2+r3+r4+r5+r6+r7+r8+r9+r10+r11+r12+r13+r14+r15+r16+r17+r18;
	u(j)=(r1*PL(j)+r2*PL(j)+r3*PL(j)+r4*PL(j)+r5*PL(j)+r6*PB(j)+r7*PB(j)+r8*PM(j)+r9*PM(j)+r10*PM(j)+r11*Z(j)+r12*PS(j)+r13*Z(j)+r14*Z(j)+r15*Z(j)+r16*NS(j)+r17*NS(j)+r18*Z(j))/rsum;
	if u(j) <= 0
       u(j) = 0;
    end
	%parameter for RFNNI
	x1 = g(j);           %input point
	x2 = u(j);
	f1(j) = exp(-1*((x1-m1(j))^2)/(2*(sigma(j)^2)) ); 
	f2(j) = exp(-1*((x1-m2(j))^2)/(2*(sigma(j)^2)) );
	f3(j) = exp(-1*((x1-m3(j))^2)/(2*(sigma(j)^2)) );
	g1(j) =  exp(-1*((x2-m4(j))^2)/(2*(sigma(j)^2)) ) ;
	g2(j) =  exp(-1*((x2-m5(j))^2)/(2*(sigma(j)^2)) ) ;
	g3(j) =  exp(-1*((x2-m6(j))^2)/(2*(sigma(j)^2)) ) ;
	fs=f1(j)+f2(j)+f3(j);
	gs=g1(j)+g2(j)+g3(j);
	%input1*input2
	mux1 = f1(j)*g1(j);   mux2 = f1(j)*g2(j);	mux3 = f1(j)*g3(j);
	mux4 = f2(j)*g1(j);   mux5 = f2(j)*g2(j);	mux6 = f2(j)*g3(j);
	mux7 = f3(j)*g1(j);   mux8 = f3(j)*g2(j);	mux9 = f3(j)*g3(j);
	mux = [mux1, mux2, mux3, mux4, mux5, mux6, mux7, mux8, mux9];
    s =sum(mux);% mux1 + mux2 + mux3+mux4 + mux5 + mux6 + mux7 + mux8 + mux9;
	%rule
	rule1(j) = w1(j)*mux1;  
    rule2(j) = w2(j)*mux2;
    rule3(j) = w3(j)*mux3;
    rule4(j) = w4(j)*mux4;
    rule5(j) = w5(j)*mux5;
    rule6(j) = w6(j)*mux6;
    rule7(j) = w7(j)*mux7;
    rule8(j) = w8(j)*mux8;
    rule9(j) = w9(j)*mux9;
	%Adjust weight
	diffjg =-(gb-g(j));
    
	%parameter for diffgu(j)
	diffg1x2 = g1(j)*(-(x2-m4(j))/sigma(j)^2);
	diffg2x2 = g2(j)*(-(x2-m5(j))/sigma(j)^2);
	diffg3x2 = g3(j)*(-(x2-m6(j))/sigma(j)^2);
	diffgx2sum = diffg1x2 + diffg2x2 + diffg3x2;
	%differential for rule1+...+rule9
	diffr147sum = (f1(j)*w1(j)+f2(j)*w4(j)+f3(j)*w7(j))*(s*diffg1x2-g1(j)*diffgx2sum*fs)/s^2;
	diffr258sum = (f1(j)*w2(j)+f2(j)*w5(j)+f3(j)*w8(j))*(s*diffg2x2-g2(j)*diffgx2sum*fs)/s^2;
	diffr369sum = (f1(j)*w3(j)+f2(j)*w6(j)+f3(j)*w9(j))*(s*diffg3x2-g3(j)*diffgx2sum*fs)/s^2;
	
	diffgu(j) = diffr147sum+ diffr258sum+ diffr369sum;
	
	%parameter for diffuw                         (question!)
	diffns=(r16+r17)/rsum;
	diffz=(r12+r14+r15+r18)/rsum;
	diffps=(r11+r13)/rsum;
	diffpm=(r8+r9+r10)/rsum;
	diffpb=(r6+r7)/rsum;
	diffpl=(r1+r2+r3+r4+r5)/rsum;
	
	diffw1 = diffns*diffgu(j)*diffjg;
	diffw2 = diffz*diffgu(j)*diffjg;
	diffw3 = diffps*diffgu(j)*diffjg;
	diffw4 = diffpm*diffgu(j)*diffjg;
	diffw5 = diffpb*diffgu(j)*diffjg;
	diffw6 = diffpl*diffgu(j)*diffjg;
	
	NS(j+1) = NS(j)+pace*(-diffw1);
	Z(j+1) = Z(j)+pace*(-diffw2);
	PS(j+1) = PS(j)+pace*(-diffw3);
	PM(j+1) = PM(j)+pace*(-diffw4);
	PB(j+1) = PB(j)+pace*(-diffw5);
	PL(j+1) = PL(j)+pace*(-diffw6);
	
    %rungekutta
    k11 = -p1*(g(j)-gb)-x(j)*g(j)+d(j);
    k21 = -p2*x(j)+p3*(i(j)-ib);
    k31 = -n*(i(j)-ib)+u(j);
    k12 = -p1*(g(j)+step/2*k11-gb)-(x(j)+step/2*k21)*(g(j)+step/2*k11)+d(j);
    k22 = -p2*(x(j)+step/2*k21)+p3*(i(j)+step/2*k31-ib);
    k32 = -n*(i(j)+step/2*k31-ib)+u(j);
    k13 = -p1*(g(j)+step/2*k12-gb)-(x(j)+step/2*k22)*(g(j)+step/2*k12)+d(j);
    k23 = -p2*(x(j)+step/2*k22)+p3*(i(j)+step/2*k32-ib);
    k33 = -n*(i(j)+step/2*k32-ib)+u(j);
    k14 = -p1*(g(j)+step*k13-gb)-(x(j)+step*k23)*(g(j)+step*k13)+d(j);
    k24 = -p2*(x(j)+step*k23)+p3*(i(j)+step*k33-ib);
    k34 = -n*(i(j)+step*k33-ib)+u(j);
    
	g(j+1) =  g(j)+step*(k11+2*k12+2*k13+k14)/6;
	x(j+1) = x(j)+step*(k21+2*k22+2*k23+k24)/6;
	i(j+1) =  i(j)+step*(k31+2*k32+2*k33+k34)/6;
	
    %start  tracking
	%Pdm for center point
	Pdm1=(x1-m1(j))/(sigma(j)^2)*exp(-((x1-m1(j))^2)/(2*(sigma(j)^2)));
	Pdm2=(x1-m2(j))/(sigma(j)^2)*exp(-((x1-m2(j))^2)/(2*(sigma(j)^2)));
	Pdm3=(x1-m3(j))/(sigma(j)^2)*exp(-((x1-m3(j))^2)/(2*(sigma(j)^2)));
	Pdm4=(x2-m4(j))/(sigma(j)^2)*exp(-((x2-m4(j))^2)/(2*(sigma(j)^2)));
	Pdm5=(x2-m5(j))/(sigma(j)^2)*exp(-((x2-m5(j))^2)/(2*(sigma(j)^2)));
	Pdm6=(x2-m6(j))/(sigma(j)^2)*exp(-((x2-m6(j))^2)/(2*(sigma(j)^2)));
	%Pdm for wdight
	Pdma1=mux1*(-(x1-m1(j))^2/sigma(j)^3-(x2-m4(j))^2/sigma(j)^3);
	Pdma2=mux2*(-(x1-m1(j))^2/sigma(j)^3-(x2-m5(j))^2/sigma(j)^3);
	Pdma3=mux3*(-(x1-m1(j))^2/sigma(j)^3-(x2-m6(j))^2/sigma(j)^3);
	Pdma4=mux4*(-(x1-m2(j))^2/sigma(j)^3-(x2-m4(j))^2/sigma(j)^3);
	Pdma5=mux5*(-(x1-m2(j))^2/sigma(j)^3-(x2-m5(j))^2/sigma(j)^3);
	Pdma6=mux6*(-(x1-m2(j))^2/sigma(j)^3-(x2-m6(j))^2/sigma(j)^3);
	Pdma7=mux7*(-(x1-m3(j))^2/sigma(j)^3-(x2-m4(j))^2/sigma(j)^3);
	Pdma8=mux8*(-(x1-m3(j))^2/sigma(j)^3-(x2-m5(j))^2/sigma(j)^3);
	Pdma9=mux9*(-(x1-m3(j))^2/sigma(j)^3-(x2-m6(j))^2/sigma(j)^3);
	Pdms=Pdma1+Pdma2+Pdma3+Pdma4+Pdma5+Pdma6+Pdma7+Pdma8+Pdma9;
	Pdmwa1=w1(j)*(Pdma1*s-mux1*Pdms)/s^2;
	Pdmwa2=w2(j)*(Pdma2*s-mux2*Pdms)/s^2;
	Pdmwa3=w3(j)*(Pdma3*s-mux3*Pdms)/s^2;
	Pdmwa4=w4(j)*(Pdma4*s-mux4*Pdms)/s^2;
	Pdmwa5=w5(j)*(Pdma5*s-mux5*Pdms)/s^2;
	Pdmwa6=w6(j)*(Pdma6*s-mux6*Pdms)/s^2;
	Pdmwa7=w7(j)*(Pdma7*s-mux7*Pdms)/s^2;
	Pdmwa8=w8(j)*(Pdma8*s-mux8*Pdms)/s^2;
	Pdmwa9=w9(j)*(Pdma9*s-mux9*Pdms)/s^2;
	Pdmy=Pdmwa1+Pdmwa2+Pdmwa3+Pdmwa4+Pdmwa5+Pdmwa6+Pdmwa7+Pdmwa8+Pdmwa9;
	wm=w1(j)*mux1+w2(j)*mux2+w3(j)*mux3+w4(j)*mux4+w5(j)*mux5+w6(j)*mux6+w7(j)*mux7+w8(j)*mux8+w9(j)*mux9;
	%y track ym
	y(j) = (rule1(j)+rule2(j)+rule3(j)+rule4(j)+rule5(j)+rule6(j)+rule7(j)+rule8(j)+rule9(j))/s;
    ym(j) = g(j);
	%Adjust width
	sigma(j+1) = sigma(j)-ata*(y(j)-ym(j))*Pdmy;
	%adjust mid
	m1(j+1) = m1(j)-ata*(y(j)-ym(j))*( (w1(j)*g1(j)+w2(j)*g2(j)+w3(j)*g3(j))*Pdm1/s-wm*gs*Pdm1/(s^2) );
	m2(j+1) = m2(j)-ata*(y(j)-ym(j))*( (w4(j)*g1(j)+w5(j)*g2(j)+w6(j)*g3(j))*Pdm2/s-wm*gs*Pdm2/(s^2) );
	m3(j+1) = m3(j)-ata*(y(j)-ym(j))*( (w7(j)*g1(j)+w8(j)*g2(j)+w9(j)*g3(j))*Pdm3/s-wm*gs*Pdm3/(s^2) );
	m4(j+1) = m4(j)-ata*(y(j)-ym(j))*( (w1(j)*f1(j)+w4(j)*f2(j)+w7(j)*f3(j))*Pdm4/s-wm*fs*Pdm4/(s^2) );
	m5(j+1) = m5(j)-ata*(y(j)-ym(j))*( (w2(j)*f1(j)+w5(j)*f2(j)+w8(j)*f3(j))*Pdm5/s-wm*fs*Pdm5/(s^2) );
	m6(j+1) = m6(j)-ata*(y(j)-ym(j))*( (w3(j)*f1(j)+w6(j)*f2(j)+w9(j)*f3(j))*Pdm6/s-wm*fs*Pdm6/(s^2) );
    %adjust weight
	w1(j+1) = w1(j)  -ata*(y(j)-ym(j))*mux1/s;
    w2(j+1) = w2(j)  -ata*(y(j)-ym(j))*mux2/s;
    w3(j+1) = w3(j)  -ata*(y(j)-ym(j))*mux3/s;
    w4(j+1) = w4(j)  -ata*(y(j)-ym(j))*mux4/s;
    w5(j+1) = w5(j)  -ata*(y(j)-ym(j))*mux5/s;
    w6(j+1) = w6(j)  -ata*(y(j)-ym(j))*mux6/s;
    w7(j+1) = w7(j)  -ata*(y(j)-ym(j))*mux7/s;
    w8(j+1) = w8(j)  -ata*(y(j)-ym(j))*mux8/s;
    w9(j+1) = w9(j)  -ata*(y(j)-ym(j))*mux9/s;
end

%print figure
figure(1)
plot( t,g(1:length(t)),'red' ); hold on
plot( t,y(1:length(t)),'blue' ); hold on
title('Glucose tracking');
xlabel('time');
legend('glucose','systemtracking');
axis([-inf, inf, 40, 100]);
figure(2)
plot( t,i(1:length(t)),'red' );hold on
title('insulin tracking');
xlabel('time');
axis([-inf, inf, 4, 14]);

figure(3)
hold on;
plot( t,u(1:length(t)),'red' );hold on
title('Control Signal');
xlabel('time');
hold off;
