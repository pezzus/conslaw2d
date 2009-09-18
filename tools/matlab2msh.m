function [] = matlab2msh(p,e,t,filename)
    %tri = [t(1:3,:)-1; t(4,:)];
    %ed = [e(1:2,:)-1; e(5,:)];
    tri = t;
    ed = e([1,2,5],:);
    pt = p;
    fid = fopen(filename,'wt');
    fprintf(fid,'# DATA\n');
    dim = [max(size(pt)),max(size(tri)),max(size(ed))];
    fprintf(fid,'%d\t%d\t%d\n',dim);
    fprintf(fid,'# POINTS\n');
    fprintf(fid,'%f\t%f\n',pt);
    fprintf(fid,'# ELEMENTS\n');
    fprintf(fid,'3\t%d\t%d\t%d\t%d\n',tri-1);
    fprintf(fid,'# EDGES\n');
    fprintf(fid,'%d\t%d\t%d\n',ed-1);
    fclose(fid);
end