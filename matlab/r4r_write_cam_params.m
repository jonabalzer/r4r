function [] = r4r_write_cam_params( filename, cam )

fid = fopen(filename, 'w');

% intrinsic params
fprintf(fid, '# dims\n');
fprintf(fid, '%d %d\n',cam.sizes);
fprintf(fid, '# focal lengths\n');
fprintf(fid, '%f %f\n',cam.f);
fprintf(fid, '# principle point\n');
fprintf(fid, '%f %f\n',cam.c);
fprintf(fid, '# radial distortion coefficients\n');
fprintf(fid, '%f %f %f %f %f\n',cam.k);
fprintf(fid, '# skew coefficient\n');
fprintf(fid, '%f\n',cam.alpha);

% extrinsic params
fprintf(fid, '# frame world -> cam \n');
fprintf(fid, '%f %f %f %f\n',cam.F(1,1),cam.F(1,2),cam.F(1,3),cam.F(1,4));
fprintf(fid, '%f %f %f %f\n',cam.F(2,1),cam.F(2,2),cam.F(2,3),cam.F(2,4));
fprintf(fid, '%f %f %f %f\n',cam.F(3,1),cam.F(3,2),cam.F(3,3),cam.F(3,4));

fclose(fid);
 
end
