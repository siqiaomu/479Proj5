
for i = 1:length(mu)
    u_mode(i,:) = Phi(:,i)*y0(i);
end

xx = xx(:,1);
mm = linspace(1,length(mu),length(mu));
[xx,mm] = meshgrid(xx,mm);


surf(xx,mm,abs(Phi'));
zlabel('\Psi','FontSize', 20);
ylabel('Mode','fontsize',18);
xlabel('x','FontSize', 18);