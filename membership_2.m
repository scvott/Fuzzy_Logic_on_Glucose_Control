%parameter of membership function for glucose
	gb = 70;a=-5; b=5; step=0.1; t=a:step:b;
	for n = 1:length(t)
	%parameter of membership function for glucose deviation
	N = [1,(-t(n)/1)]; ZERO = [(t(n)+1),(1-t(n))];P = [t(n)/1,1];
	%membership function for glucose deviation
	neg(n) = max(min(N),0);
	zero(n) = max(min(ZERO),0);
	pos(n) = max(min(P),0);
	end
	figure(1)
    %title("Glucose Deviation Membership Function");
	xlabel(' error differential'); 
    ylabel('degree of membership');
    hold on
	plot(t,neg(1:length(t)));
	plot(t,zero(1:length(t)));
	plot(t,pos(1:length(t)));
	hold off
    legend('Negative','Zeor','Positive'); 