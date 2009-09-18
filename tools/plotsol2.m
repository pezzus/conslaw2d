function [] = plotsol2(directory,fileavi)
    N = 100;
    for i = 1:N
        filename = [directory, 'solution', num2str(i-1,'%0.4d'), '.dat'];
        data = load(filename);
        % Numero di variabili
        x = data(:,1);
        y = data(:,2);
        z = data(:,3);
        [XI,YI] = meshgrid(min(x):0.02:max(x), min(y):0.02:max(y));
        ZI = griddata(x,y,z,XI,YI);
        %trisurf(t(1:3,:)',p(1,:),p(2,:),z);
        shading interp;
        colormap cool;
        %vertices = [x, y, z];
        %connect = t(1:3)';
        %patch('Vertices',vertices,'Faces',connect);
        %contourf(XI,YI,ZI1,10,'LineStyle','none');
        %surf(XI,YI,ZI1,'LineStyle','none');
        surfl(XI,YI,ZI);
        shading interp
        colormap(gray);
        %axis([-1,1,-1,1,0,10]);
        %caxis([0,10]);
        %view(-0.5,90);
        M(i) = getframe;
        %hold on;
        %quiver(XI,YI,ZI2,ZI3); axis image; hold off;
        %pause;
    end
    %movie(M);
    %movie2avi(M,fileavi);
end