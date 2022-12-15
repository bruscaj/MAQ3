function [NMatSDD,CMatSDD,Points,Interior,Boundary,weight,Dvvs] = slabBdryCubeDisc(depth,x0,x1,N)
% depth = depth of stencil
% [x0,x1] = cube dimension *Note line 115*
% N = # of grid points

% Note: This code assume an Linf stencil and makes use of specific
% structure related to this stencil. Building a more general code would
% either knowing similar properties for each possible stencil type, or
% would require finding quadrature points via solving a system of equations

% Initialize stencil array
% depth = 5; % depth of the stencil
D = vector_creation(depth,'inf',3); % creates points in an Linf Stencil
DD = [D D]; % Create structure to reorganize stencil points 
DD(:,1:2:end) = D; % Original stencil points
DD(:,2:2:end) = -D; % Layer in opposite stencil points to be adjacent
% DD will have ~double of each stencil point
DD = unique(DD','rows','stable'); % Drop repeats, keep structure
D = DD'; % Reassign refomatted stencil array to D
dCount = length(D); % Number of stencil points

% Calculate Quadrature weights
F1 = D(1,:) == depth; % Stencil Points on the first face, x = depth
F2 = D(2,:) == depth; % Stencil Points on the second face, y = depth
F3 = D(3,:) == depth; % Stencil Points on the third face, z = depth
Fcs = [F1;F2;F3]; % Store as a single array for easier access
Wcs = zeros(size(Fcs)); % Will calculate weights on each face
% Our ultimate goal is to find the quadrature weight associate to each pair
% of stencil points that form a SDD. We will find the weight on each face,
% then sum those weights together. We will treat the quadrature as a single
% sum over all points **(This may be slower)**

% Loop through each face
for kk = 1:3
    dd = D(setdiff(1:3,kk),Fcs(kk,:))/depth; % Normalized stencils points 
    gp = sqrt(length(dd)); % number of grid points **(assumes Linf)**
    s = linspace(-1,1,gp); t = s; % redisc for structure
    % Rather than trying to use inherent strucutre, we will create an
    % identical array of stencil points which follows an obvious structure,
    % then map the quadrature weights here onto the original array
    ws = quadWeights(s,2); wt = ws; % simp. weights in each direction 
    wst = ws*wt'; % simp. over face by outter product
    [S,T] = meshgrid(s,t); % generate grid for structure
    ST = 1./sqrt(1+S.^2+T.^2).^3; % Calculate weight function
    W = wst.*ST; % array of weights from simp. and cube face weight
    d = [S(:)';T(:)']; % grid array with imposed structure
    w = W(:)'; % weight array with impose structure
    [~,indkk] = ismember(round(dd'*depth),round(depth*d'),'rows'); 
    % find mapping from generate structure back to original
    % this works because d(:,indkk) = dd(:)
    w = w(indkk); % map weights to correct nodes
    Wcs(kk,Fcs(kk,:)) = w; % assign correct weights to corresponding face
end
weight = (sum(reshape(sum(Wcs),2,[])))'; % sum weights from all faces, reshape
% to have opposite nodes together, sum over opposite nodes
% weight 

% Things to note:
% 1) Would probably be faster to use some inherent structure rather than
% ismember. This can be improved at some point, but for now this is fine
% 2) The weights should sum to 2pi, but we're getting something slighly
% different. The error seems to be two orders of magnitude less than the
% individual weights, so it doesn't seem likely that something is being
% added wrong. On the otherhand, if it were just a rounding error, it seems
% more likely that the error would be of a much smaller scale. Will need to
% figure out what is causing this issue
% 3) The way we are storing the weights may lose some symmetry which could
% potentially be exploited for faster computations in the future. 
% 4) Loop has some redundancies, and is slower than necessary. This should
% be fixed whenever we use the inherent structure to compute weights faster

% Thoughts: for 1)&4), could probably just create the D array while we
% create the quadarture weights, though it may take some work to do this
% efficiently. Fine to use this form for now.
% Thoughts: for 2) this may be from the way we are applying the weight
% function. Specifically, the error may come from the quadrature being
% applies to the weight function. 

% Note: we are going to build a discretization array which we can use for
% calculating finite difference on a cartesian grid. For now, we will
% assume a cube domain with equal spacing in each direction. In order to
% impose the wide stencils near the boundary, we will assume that the
% solution data is know at and *near* the boundary. Futute implementations
% will require some additional method to oversample the boundary.

% Cube discretization 
% N = 2^4+1; % points along each direction
% x0 = 0; x1 = 1; % x boundaries
% y0 = 0; y1 = 1; % y boundaries
% z0 = 0; z1 = 1; % z boundaries

dtol = depth+1; % distance to overshoot boundary for wide stencils
h = (x1-x0)/(N-1); % grid resolution, assuming uniform in each direction
x= (x0-dtol*h):h:(x1+dtol*h); % h res. sampling of [x0-dtol*h,x1+dtol*h]
% y = x; z = x; % assuming x = y = z 
[xx,yy,zz] = ndgrid(x); % create grid points
X = xx(:); Y = yy(:); Z = zz(:); % map to linear index/point cloud
Ind_Temp = (1:length(X))'; % temporary index of grid points
% The stencil D is in terms of the triple index(i,j,k), but we want to
% apply to the linear index l

% Creating Neighbor matrix
ind = 1:length(x); % sampling index of x,y,z
[ix,iy,iz] = ndgrid(ind); % grid of index values (i,j,k)
Ix = ix(:); Iy = iy(:); Iz = iz(:); % squish index triples to index cloud
I1 = ones(size(Ix)); % array of 1s for regression
II = [I1 Ix Iy Iz]; % organize values into single array
map = II\Ind_Temp; % regression to find map from indices to cloud index
map = round(map); % map will have integer coefficients
% We know there is a map such that L(i,j,k) = a*i+b*j+c*k+d which comes
% from how we squish the data. Could find it using the inherent structure
% of the squish we use, but this fine for a first implementation

Points_Old = [X Y Z]; % Original Points cloud
Ext_Old = find(X >x1 | X <x0 | Y > x1 | Y < x0 | Z > x1 | Z <x0); 
% Index of exterior points
Int_Temp = setdiff(Ind_Temp,Ext_Old); % Points not in the exterior
NMat = repmat(Int_Temp,1,dCount); % initialize neighbor matrix for non ext
% NMat has the structure that each row corresponds to a point in the grid,
% and each column corresponds to a direction in the stencil D. The value is
% the index of the neighbor.
% **(NMat will initital contain points outside of Ind_Temp)**
dL = (D'*map(2:end))'; % Convert stencil difference to index cloud diff
NMat = NMat+repmat(dL,length(Int_Temp),1); % Assign values to NMat

vCount = length(weight); % Number of unique directions
NMat_Temp = zeros(length(Int_Temp),(vCount)*3); % Initialize SDD structure 
% the SDD structure groups each 3 points in the centered diff
NMat_Temp(:,1:3:(vCount)*3) = repmat(Int_Temp,1,length(1:3:(vCount)*3)); 
% Assign self point
NMat_Temp(:,2:3:(vCount)*3) = NMat(:,1:2:end); % assign 'forward' point
NMat_Temp(:,3:3:(vCount)*3) = NMat(:,2:2:end); % assign 'backward' point
% NMat_Temp has neighbors which are in the exterior, but we only want
% neighbors up to the boundary

[tempind,~] = find(ismember(NMat_Temp,Ext_Old)); % Find ext neighbor row
%These points will make up our boundary
Bdry_Old = NMat_Temp(unique(tempind),1); % Original index of boundary nodes
Int_Old = setdiff(Int_Temp,Bdry_Old); % Original index of interior nodes
NMat_Old = NMat_Temp(ismember(NMat_Temp(:,1),Int_Old),:); % Restrict to int
% Would like to map original indices to a new point cloud without exterior
% points. Will impose structure such that Indices take the form [Int;Bdry]

intLen = length(Int_Old); bdryLen = length(Bdry_Old); % Number of cor. pts
Interior = (1:intLen)'; Boundary = (intLen+(1:bdryLen))'; % Generate structure
old2new = zeros(size(Ind_Temp)); % Inititalize map from original to new pts
old2new(Int_Old) = Interior; % Maps old interior index to new int index
old2new(Bdry_Old) = Boundary; % Map old bdry index to new bdry index
NMatSDD = old2new(NMat_Old); % Apply map to NMat_Old
Points= [Points_Old(Int_Old,:); Points_Old(Bdry_Old,:)]; % New pt cloud

% Coefficient Matrix
CMatSDD = zeros(size(NMatSDD)); % Inititalize CMatSDD
Coef = vecnorm(h*D).^-2; % Vector of SDDs coef
% This relies on the fact that we are not oversampling the boundary yet,
% but instead are using exact centered difference everywhere
CMatSDD(:,2:3:end) = repmat(Coef(1:2:end),length(Interior),1); % Frwd coef
CMatSDD(:,3:3:end) = repmat(Coef(2:2:end),length(Interior),1); % Bkwd coef
CMatSDD(:,1:3:end) = -(CMatSDD(:,2:3:end)+CMatSDD(:,3:3:end)); % Cent coef

% SDD Matrices
Dvvs = cell(vCount,1); % Init cell of SDD mats
% Loop through each direction
for i = 1:vCount
    % Builds each sparse matrix
   Dvvs{i} = sparse( repmat(Interior,1,3), [NMatSDD(:,i*3-2) NMatSDD(:,i*3-1) NMatSDD(:,i*3)], [CMatSDD(:,i*3-2) CMatSDD(:,i*3-1) CMatSDD(:,i*3)], length(Interior), length(Points));
end

end