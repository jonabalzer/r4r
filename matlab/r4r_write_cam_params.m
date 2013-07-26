function [] = r4r_write_cam_params( filename, cam )

fileOut = fopen(filename, 'w');

% intrinsic params
fprintf(fileOut, '# dims\n');
fprintf(fileOut, '%d %d\n',cam.sizes);
fprintf(fileOut, '# focal lengths\n');
fprintf(fileOut, '%f %f\n',cam.f);
fprintf(fileOut, '# principle point\n');
fprintf(fileOut, '%f %f\n',cam.c);
fprintf(fileOut, '# radial distortion coefficients\n');
fprintf(fileOut, '%f %f %f %f %f\n',cam.k);
fprintf(fileOut, '# skew coefficient\n');
fprintf(fileOut, '%f\n',cam.alpha);

% extrinsic params
fprintf(fileOut, '# frame world -> cam \n');
fprintf(fileOut, '%f %f %f %f\n',cam.F(1,1),cam.F(1,2),cam.F(1,3),cam.F(1,4));
fprintf(fileOut, '%f %f %f %f\n',cam.F(2,1),cam.F(2,2),cam.F(2,3),cam.F(2,4));
fprintf(fileOut, '%f %f %f %f\n',cam.F(3,1),cam.F(3,2),cam.F(3,3),cam.F(3,4));

fclose(fileOut);
 
end
