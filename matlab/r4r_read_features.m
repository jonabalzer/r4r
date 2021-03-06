function [ obj ] = r4r_read_features( filename )


fid = fopen(filename, 'r');

if(fid<0)  
   error('ERROR: file not found!');    
end

% read comment
comment = fgets(fid);

% read no of features
n = fscanf(fid, '%d\n',1);
features = cell(1,n);

for k=1:n

feature = struct;    
    
% read pixel position
u = fscanf(fid, '%f',1);
v = fscanf(fid, '%f\n',1);
feature.location = [u,v];

% read scale
feature.scale = fscanf(fid, '%d\n',1);

% read quality
feature.quality = fscanf(fid, '%f\n',1);
 
% number of descriptors
nod = fscanf(fid, '%d\n',1);
descriptors = cell(1,nod);

for j=1:nod

     descriptor = struct;
    
     % name and size
     descriptor.name = fscanf(fid, '%s\n',1);
     
     % data
     descriptor.data = r4r_read_dense_mat_from_stream(fid);
         
     descriptors{j} = descriptor;     
     
end

feature.descriptors = descriptors;

features{k} = feature;

end

obj.comment = comment;
obj.features = features;


fclose(fid);

end

