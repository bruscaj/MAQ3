%Matt Cassini
%Creating vectors of identical stencil width for quadrature over 2D or 3D
%Prof Hamfeldt

%Creates vectors of given dimension n = 2,3, stencil width w, and norm type
%integer 1 or string 'inf'
%Returns nx2 matrix v sorted according to increasing theta vector t for
%2D
%Returns nx3 matrix v sorted according to increasing theta vector t, vector
%phi, and magnitude rho using 2-norm for 3D

%Comment out inputs if running as function

function v = vector_creation(w,norm,n)
%clear; clc;
%w = 3;
%n = 3;
%norm = 'inf'; norm = 1;

p = 0; %Initializing phi,rho for 3d
r = 0;
if n == 2 %2D vector creation
    
    if norm == 1
        %width = i + j
        v = [];
        [x,y] = meshgrid(-w:w,-w:w); %cartesian product of all possible elements for range 0 to pi for theta
        cartprodxy = [x(:) y(:)];
        v(1:2,:) = transpose(cartprodxy);
        s = sum(abs(v));
        v = v(:,s==w); %deletes all vectors not meeting criteria (summing to w)
        
    elseif norm == 'inf'
        %width = max(i,j)
        v = [];
        [x y] = meshgrid(-w:w,-w:w);
        cartprodxy = [x(:) y(:)];
        v(1:2,:) = transpose(cartprodxy);
        m = max(abs(v));
        v = v(:,m==w); %deletes all vectors not meeting criteria (max = w)
        
        v(:,v(1,:)==w & v(2,:)==0) = []; %remove endpoint, not necessary for periodic function
    else
        disp('Norm not supported'); %For our purposes, only want 1-norm and infinity-norm
    end
    v = transpose(unique(transpose(v),'rows')); %deletes duplicates
elseif n == 3 %3D vector creation
    if norm == 1
        %width = i + j
        v = [];
        [x,y,z] = ndgrid(-w:w,-w:w,-w:w); %cartesian product of all possible elements for range 0 to pi for theta and phi
        cartprodxyz = [x(:) y(:) z(:)];
        v(1:3,:) = transpose(cartprodxyz);
        s = sum(abs(v));
        v = v(:,s==w); %deletes all vectors not meeting criteria (max = w)
         
        
    elseif norm == 'inf'
        %width = max of i,j,k
        v = [];
        [x,y,z] = ndgrid(-w:w,-w:w,-w:w); %cartesian product of all possible elements for range 0 to pi for theta and phi
        cartprodxyz = [x(:) y(:) z(:)];
        v(1:3,:) = transpose(cartprodxyz);
        m = max(abs(v));
        v = v(:,m==w); %deletes all vectors not meeting criteria (max = w)
        
        %remove duplicates (negative and positive of vector acts the same)
        ind1 = v(1,:)==-w & v(2,:)==0 & v(3,:)==0;
        v(:,ind1) = [];

        ind2 = v(1,:)==0 & v(2,:)==-w & v(3,:)==0;
        v(:,ind2) = [];

        ind3 = v(1,:)==0 & v(2,:)==0 & v(3,:)==-w;
        v(:,ind3) = [];
        
        %v(v(1,:)==w & v(2,:)==0 &v(3,:)==0) = []; %remove endpoint, not necessary for periodic function
    else
        disp('Norm not supported'); %For our purposes, only want 1-norm and infinity-norm
    end
    v = transpose(unique(transpose(v),'rows'));
else
    disp('Dimensions not supported'); %Only going up to 3D
end

end