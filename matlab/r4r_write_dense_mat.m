function [ ] = r4r_write_dense_mat( filename, A )


fid = fopen(filename, 'w');

if(fid<0)  
   error('ERROR: file not found!');    
end



if isa(A,'double')
    type = 11;
    data = A(:);
end

if isa(A,'double') && size(A,3)==3
    type = 111;

    data = zeros(1,numel(A));
    x = A(:,:,1);
    y = A(:,:,2);
    z = A(:,:,3);
    data(1:3:length(data)) = x(:);
    data(2:3:length(data)) = y(:);
    data(3:3:length(data)) = z(:);

end

fprintf(fid, '%d %d %d\n',size(A,1),size(A,2),type);


fwrite(fid,data,'double');

fclose(fid);

end

