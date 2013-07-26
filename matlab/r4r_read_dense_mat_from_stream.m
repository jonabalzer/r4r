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

if(type==1)
    data = fread(fid,[nrows,ncols],'uint8');
end


end

