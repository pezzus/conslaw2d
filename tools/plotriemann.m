N = 100;
Nt = 10;
x = linspace(-1,1,N);
tt = linspace(0,1,Nt);
for j = 1:Nt
    QQ = [];
    for i = 1:N
        Q = RiemannProblem(x(i),tt(j),4,0,1.6,1,0,0.4);
        QQ = [QQ Q'];
    end
    plot(x,QQ);
    pause;
end