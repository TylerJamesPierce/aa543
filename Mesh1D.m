function [index,x,dx] = Mesh1D( range, imax, xc, SF )
%MESH1D  create a 1xN continuous non-uniform grid or mesh

%Tyler James Pierce
%tjp644@u.washington.edu

%Version History    01/14/14: Created
%                   01/20/14: Added non-uniform hyperbolic tangent
%                   01/22/14: Added non-uniform polynomial

%   Inputs:
%   range - spatial range for verticies to be located within
%   imax - number of verticies to be generated
%   xc - the location to cluster the grid around
%   SF - The scaling factor (dx_max/dx_min) Typically [1.1-2]

ni = nargin; %check number of inputs
n=imax;
if ni<2
    error('invalid number of inputs arguments');
end
xmin=min(range);
xmax=max(range);
L=xmax-xmin; %total length
dxconstant=L/(n-1);
xavg=xmin:dxconstant:xmax;

if SF==1 % uniform grid (dx=constant
    index=1:n;
    dx=L/(n-1)*ones(1,n);
    x=xmin:L/(n-1):xmax;
elseif SF>1 %grid clustering
    index=1:n;
    SF=dxconstant/SF;
    syms ic %define a symbolic ic to solve for the integer ic of xc.
    b=solve((xc-xmin-ic*SF)/(L-n*SF)==(2*ic^3-5*ic^2+3*ic)/(2*n^3+...
            3*n^2-6*ic*n*(1+n)+6*ic^2*n));
    % check all the solutions to the cubic equation and take the real value
    for i=1:numel(b)
        if isreal(b(i))==1
            iC=double(b(i)); % Define the psuedocluster index
            break
        end
    end
    %Solve for the slope dx=k*(i-ic)^2+SF
    k=(L-(n)*SF)/(iC^2*n-iC*n*(1+n)+(n*(n+1)*(2*n+1))/6);
    if k<0
        error('negative dx, check inputs');
    end

    dx=ones(1,n-1);x=ones(1,n); %initialize vectors
    x(1)=xmin;
    for int=1:n-1
        dx(int)=n/(n-1)*(k*(int-iC)^2+SF);
        x(int+1)=x(int)+dx(int);
    end
end
end
% Old versions
%
% This is the hyperbolic tangent function
% This is only piecewise continuous. It is continuous on both sides of
% xcluster, more fitted for progression outward from a body.
% if hyperbolic tangent parameter is turned on:
% xlow=xavg(1:iSource);
% xhigh=xavg(iSource+1:xmax);
% nHigh=nVerticies-iSource-1;%number of indexes from xsource to xmax
% nLow=iSource;
% sLow=normtanh(SF,nLow);
% sHigh=normtanh(SF,nHigh+2);
% sNew=sHigh(2:end);              
% x=zeros(1,nVerticies);
% x(1:nLow)=wrev(xlowr);
% xUppr=xc.*ones(1,nHigh+1)+sNew(1:nHigh+1).*(xmax-xc);
% x(nLow+1:nVerticies)=xUppr;
% dx=ones(1,nVerticies-1);
% for i=1:nVerticies-1
%     dx(i)=x(i+1)-x(i);
% end
% % End of the hyperbolic tangent function
% 
% % This is the polynomial fitted grid clustering
% % This is completely continuous but may be steeper than desired for
% % normal grid distributions.
% %   (xc-xmin-ic*SF)/(xmax-xmin-N*SF)=
% %   (2*ic^3-5*ic^2+3*ic)/(2*N^3+3*N^2-6*ic*(N+N^2)+6*ic^2*N)
% %   
% %   dxmax/dxmin = { k/SF*ic^2+1         if xc<L/2
% %                   k/SF*(N-ic)^2+1     if xc>L/2
% %
% dx(i)=K(i-ic)^2+SF; %dx exists in [1,N]
% bnum=-1*(3*xc*(SF-1)*xmax);
% if xc<=L/2
%     cf=K/SF*(ic^2)+1;
%     
%     b=bnum/(n*(3-3*xc*(2+n*(SF-1))+n.^2*(SF-1)+3*SF*xc.^2));
%     c=b*((2*xc-1)*n-xc.^2*SF)/(xc*(SF-1));
% elseif xc>=L/2
%     cf=K/SF*(N-ic)^2+1;
%     
%     b=bnum/(n*(-3*n*(-1+xc+xc*SF)+n.^2*(SF-1)+3*SF*xc.^2));
%     c=b*((2*xc-1)-xc.^2*SF)/(xc*(SF-1));
% else
%     error('invalid input');
% end
% a=-b/(3*xc);
% 
% delta_x=zeros(nVerticies-1,1);
% x=zeros(nVerticies,1);
% x(1)=xmin;x(nVerticies)=xmax;
% %SF=SF*L;
% % editable parameters:
% % constantDelta_x=2;
% % if constantDelta_x==1
% %     delta_x = L/(nVerticies-1)*ones(nVerticies-1); %constant vector
% % elseif constantDelta_x==2
% %     for iter=1:nVerticies-1
% %         delta_x(iter) = L/(nVerticies-1)/(1+SF*exp(-L^2*(x(iter)-xSource)^2/(2)));%den = 2*c^2
% %         x(iter+1)=x(iter)+delta_x(iter);
% %         
% %     end
% % end
% % xerror=xmax-x(end);
% for iter=2:nVerticies-1
%     %delta_x(iter) = L/(nVerticies-1)/(1+L*SF*exp(-L^2*(x(iter)-xSource)^2/(2)))+xerror/(nVerticies);
%     x(iter)=a*iter^3+b*iter^2+c*iter;
% end
% index=1:nVerticies;
% 
% end
