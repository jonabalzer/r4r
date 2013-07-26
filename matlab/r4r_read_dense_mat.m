function [ data ] = r4r_read_dense_mat( filename )

fid = fopen(filename, 'r');

if(fid<0)  
   error('ERROR: file not found!');    
end

data = r4r_read_dense_mat_from_stream(fid);

fclose(fid);

end

