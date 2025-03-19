

function [LSP,LCP] = CMDSP(varargin) % image,radius,neighbors,mapping,mode)
  
   % Check number of input arguments.
   error(nargchk(1,7,nargin));
   image=varargin{1};
   d_image=double(image);
  
   if nargin==1
      spoints=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
      neighbors=8;
      mapping=0;
      mode='h';
   end
   if (nargin == 2) && (length(varargin{2}) == 1)
       error('Input arguments');
   end
   if (nargin > 2) && (length(varargin{2}) == 1)
      radius=varargin{2};
      neighbors=varargin{3};
      spoints=zeros(neighbors,2);
      
      a = 2*pi/neighbors;
     
      for i = 1:neighbors
          spoints(i,1) = -radius*sin((i-1)*a);
          spoints(i,2) = radius*cos((i-1)*a);
      end
      if(nargin >= 4)
        mapping=varargin{4};
          if(isstruct(mapping) && mapping.samples ~= neighbors)
             error('Incompatible mapping');
          end
      else
          mapping=0;
      end
      if(nargin >= 6)
         mode=varargin{5};%
         thre = varargin{6};
         thre_M = varargin{7};
      else
         mode='h';
      end
   end

    if (nargin > 1) && (length(varargin{2}) > 1)
    spoints=varargin{2};
    neighbors=size(spoints,1);
    
    if(nargin >= 3)
        mapping=varargin{3};
        if(isstruct(mapping) && mapping.samples ~= neighbors)
            error('Incompatible mapping');
        end
    else
        mapping=0;
    end
    
    if(nargin >= 4)
        mode=varargin{4};
    else
        mode='h';
    end   
    end


  [ysize xsize] = size(image);
  miny=min(spoints(:,1));
  maxy=max(spoints(:,1));
  minx=min(spoints(:,2));
  maxx=max(spoints(:,2));
  % Block size, each LBP code is computed within a block of size bsizey*bsizex
  bsizey=ceil(max(maxy,0))-floor(min(miny,0))+1;
  bsizex=ceil(max(maxx,0))-floor(min(minx,0))+1;
  % Coordinates of origin (0,0) in the block
  origy=1-floor(min(miny,0));
  origx=1-floor(min(minx,0));
  % Minimum allowed size for the input image depends
  % on the radius of the used LBP operator.
  if(xsize < bsizex || ysize < bsizey)
     error('Too small input image. Should be at least (2*radius+1) x (2*radius+1)');
  end
  % Calculate dx and dy;
  dx = xsize - bsizex;
  dy = ysize - bsizey;
  % Fill the center pixel matrix C.
  C = image(origy:origy+dy,origx:origx+dx);  
  d_C = double(C);  
  bins = 2^neighbors;

  LSP=zeros(dy+1,dx+1);
  LCP=zeros(dy+1,dx+1);
  %Compute the LBP code image
   for i = 1:neighbors
       y = spoints(i,1)+origy;
       x = spoints(i,2)+origx;
       % Calculate floors, ceils and rounds for the x and y.
       fy = floor(y); cy = ceil(y); ry = round(y);
       fx = floor(x); cx = ceil(x); rx = round(x);
       % Check if interpolation is needed.
       if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
       % Interpolation is not needed, use original datatypes
           N = d_image(ry:ry+dy,rx:rx+dx);
           D{i} = N >= d_C;   

           Diff{i} = 1./(1+exp(N-d_C)); 
           MeanDiff(i) = mean(mean(Diff{i}));
       else
       %Interpolation needed, use double type images 
           ty = y - fy;
           tx = x - fx;
       % Calculate the interpolation weights.
           w1 = (1 - tx) * (1 - ty);
           w2 =      tx  * (1 - ty);
           w3 = (1 - tx) *      ty ;
           w4 =      tx  *      ty ;
       % Compute interpolated pixel values
           N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
               w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
           D{i} = N >= d_C;     

           Diff{i} = 1./(1+exp((N-d_C)));
           MeanDiff(i) = mean(mean(Diff{i}));
      end  
end

DiffThreshold = mean(MeanDiff);

for i=1:neighbors
  % Update the result matrix.
  
  v = 2^(i-1);

  LSP = LSP + v*(Diff{i}>=thre_M);
end

threC = thre(origy:origy+dy,origx:origx+dx);   
LCP = d_C>=threC;

%Apply mapping if it is defined
if isstruct(mapping)
    bins = mapping.num;
     sizarray = size(LSP);

    LSP = LSP(:);

    LSP = mapping.table(LSP+1);

    LSP = reshape(LSP,sizarray);

end

if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh'))

    LSP=hist(LSP(:),0:(bins-1));
    if (strcmp(mode,'nh'))

        LSP=LSP/sum(LSP);
    end
else

end






