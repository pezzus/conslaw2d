function [] = plotsol(directory,N,cax)
    for i = 1:N
        filename = [directory, 'solution', num2str(i-1,'%0.4d'), '.dat'];
        data = load(filename);
        % Numero di variabili
        x = data(:,1);
        y = data(:,2);
        rho = data(:,3);
        velx = data(:,4);
        vely = data(:,5);
        pres = data(:,6);
        [XI,YI] = meshgrid(min(x):0.005:max(x), min(y):0.005:max(y));
        % Plot densita'
        ZI = griddata(x,y,rho,XI,YI);
        %s1 = subplot(2,2,1);
        contour(XI,YI,ZI,10);
        %caxis(cax(1,:));
        s1 = title('Densità');
        % Plot pressione
        %ZI = griddata(x,y,pres,XI,YI);
        %s2 = subplot(2,2,2);
        %contourf(XI,YI,ZI,10);
        %caxis(cax(2,:));
        %s2 = title('Pressione');
        % Plot velocita'
        %ZI = griddata(x,y,sqrt(velx.*velx+vely.*vely),XI,YI);
        %s3 = subplot(2,2,3:4);
        %contourf(XI,YI,ZI,10,'LineStyle','none'); hold on;
        %caxis(cax(3,:));
        %[DX,DY] = gradient(ZI,.5,.5);
        %quiver(XI,YI,DX,DY); hold off;
        %s3 = title('Velocità');
        %plot(XI(49,:),ZI(49,:));
        %axis([-1,1,9,12]);
        %surf(XI,YI,ZI1,'LineStyle','none');
        %surfl(XI,YI,ZI);
        %shading interp
        %colormap(gray);
        %axis([-1,1,-1,1,0.9,2.1]);
        %caxis([0.9,2.1]);
        M(i) = getframe;
        %hold on;
        %quiver(XI,YI,ZI2,ZI3); axis image; hold off;
        %pause;
    end
    %movie(M);
    %movie2avi(M,fileavi);
end