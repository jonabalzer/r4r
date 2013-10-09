function [ data ] = r4r_read_dense_mat_from_stream( fid )

nrows = fscanf(fid, '%d',1);
ncols = fscanf(fid, '%d',1);
type = fscanf(fid, '%d\n',1);

if(type==11)
    data = fread(fid,[nrows,ncols],'double');
end

if(type==8)
    data = fread(fid,[nrows,ncols],'float');
end

if(type==6)
    data = fread(fid,[nrows,ncols],'int32');
end

if(type==1 || type==2)
    data = fread(fid,[nrows,ncols],'uint8');
end


if(type==111)
    data = fread(fid,[1,nrows*ncols*3],'double');
    
    x = reshape(data(1:3:length(data)),[nrows,ncols]);
    y = reshape(data(2:3:length(data)),[nrows,ncols]);
    z = reshape(data(3:3:length(data)),[nrows,ncols]);
    
    data = cat(3,x,y,z);
end


end

