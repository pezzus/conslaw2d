data = load('../test/sodproblem/data/solution0105.dat');
% Riordino
data = sortrows(data);
%t = 0.00505445;
%t = 0.0499941;
%t = 0.198247;
t = 0.36638;
x = data(:,1);
Qapprox_rho = data(:,2)';
Qapprox_u = data(:,3)';
Qapprox_p = data(:,4)';
Qexact_rho = [];
Qexact_u = [];
Qexact_p = [];
for i = 1:size(data,1)
    Q = RiemannProblem(x(i),t,1.0,0.75,1.0,0.125,0.0,0.1);
    Qexact_rho = [Qexact_rho Q(1)];
    Qexact_u = [Qexact_u Q(2)];
    Qexact_p = [Qexact_p Q(3)];
end
subplot(3,1,1);
h1 = plot(x,Qexact_rho,'--o','LineWidth',1.0,'Color','k','MarkerFaceColor','w');
hold on;
plot(x,Qapprox_rho,'-sq','LineWidth',2.0,'Color','k','MarkerFaceColor','w');
axis([-1,1.05,0,1.1]);
set(gca,'xtick',[])
h1 = title('Densita''');
subplot(3,1,2);
h2 = plot(x,Qexact_u,'--o','LineWidth',1.0,'Color','k','MarkerFaceColor','w');
hold on;
plot(x,Qapprox_u,'-sq','LineWidth',2.0,'Color','k','MarkerFaceColor','w');
axis([-1,1.05,-0.1,2]);
set(gca,'xtick',[])
h2 = title('Velocita''');
subplot(3,1,3);
h3 = plot(x,Qexact_p,'--o','LineWidth',1.0,'Color','k','MarkerFaceColor','w');
hold on;
plot(x,Qapprox_p,'-sq','LineWidth',2.0,'Color','k','MarkerFaceColor','w');
axis([-1,1.05,0,1.1]);
set(gca,'xtick',[])
h3 = title('Pressione');
set(gca,'xtickMode', 'auto')
% Norme
dx = x(2:end)-x(1:end-1);
dQ_rho = abs(Qexact_rho(2:end) - Qapprox_rho(2:end));
dQ_u = abs(Qexact_u(2:end) - Qapprox_u(2:end));
dQ_p = abs(Qexact_p(2:end) - Qapprox_p(2:end));
NormL1_rho = sum(dx'.*dQ_rho);
NormL1_u = sum(dx'.*dQ_u);
NormL1_p = sum(dx'.*dQ_p);
