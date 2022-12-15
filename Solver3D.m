MM = 4;
E = zeros(1,MM);
H = zeros(1,MM);
for kk = 1:MM
x0 = 0; x1 = 1;
N = 2^(kk+1)+1;
h = (x1-x0)/(N-1); epsilon = h^2;
depth = floor(h^(-1/3));
[NMatSDD,CMatSDD,Points,Interior,Boundary,weight,Dvvs] = slabBdryCubeDisc(depth,x0,x1,N);
vCount = length(weight);
uExact = @(x,y,z) 1/2*(x.^2+y.^2+z.^2);
fExact = @(x,y,z) ones(size(x));


uTrue = uExact(Points(:,1),Points(:,2),Points(:,3));
uSoln = [zeros(size(Interior));uExact(Points(Boundary,1),Points(Boundary,2),Points(Boundary,3))];
F = fExact(Points(Interior,1),Points(Interior,2),Points(Interior,3));
aproxMAOp = @(u) 4*pi^2*(SDDMat(NMatSDD,CMatSDD,u,vCount,epsilon).^(-3/2)*weight).^-2+min(min(SDDMat(NMatSDD,CMatSDD,u,vCount,-Inf),epsilon),[],2);
bZ = zeros(size(Boundary));

resid = norm(aproxMAOp(uSoln)-F,'inf');
finished = true;

max_count = 100;
stepcount = 0;

uNewt = uSoln;

while resid > h^2
deltaU = [newtUpdate3D(NMatSDD,CMatSDD,Dvvs,weight,uNewt,F,epsilon);bZ];
    % Only updated in the Interior
    
    alpha = 1;
    uNewtTemp = uNewt -alpha*deltaU;
    residTemp = norm(aproxMAOp(uNewtTemp)-F,inf);
    % alpha is our dampening coeffecient, reset it to 1 at each iteration,
    % newtUpdate_modQuad takes in the current u and the values of F at the interior
    % and gives an update for newton's method. uNewtTemp and residTemp are
    % the values of u and the residual with this update
    
    while residTemp > resid
        
        alpha = alpha/2;
        uNewtTemp = uNewt -alpha*deltaU;
        residTemp = norm(aproxMAOp(uNewtTemp)-F,inf);
        % We bisect alpha until the residual decreases
        
        if alpha < 10^-16
            finished = false;
            break
        end
        
    end
    if ~finished
        break
    end
    resid = residTemp;
    uNewt = uNewtTemp;
    stepcount = stepcount + 1;
    if stepcount >= max_count
        break
    end
end
E(kk) = norm(uNewt-uTrue,inf)
H(kk) = h
end
