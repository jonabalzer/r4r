function [ cam ] = r4r_read_cam_params( filename )

fileIn = fopen(filename, 'r');

if(fileIn<0)  
   error('ERROR: file not found!');    
end

fgets(fileIn);

sizes = fscanf(fileIn, '%d',2);
 
fgets(fileIn);
fgets(fileIn);

f = fscanf(fileIn, '%f',2);

fgets(fileIn);
fgets(fileIn);

c = fscanf(fileIn, '%f',2);

fgets(fileIn);
fgets(fileIn);

k = fscanf(fileIn, '%f',5);

fgets(fileIn);
fgets(fileIn);

alpha = fscanf(fileIn, '%f',1);

fgets(fileIn);
fgets(fileIn);

F = fscanf(fileIn, '%f',[3,4]);
cam.F = F';

cam.sizes = sizes;
cam.f = f;
cam.c = c;
cam.k = k;
cam.alpha = alpha;

end

