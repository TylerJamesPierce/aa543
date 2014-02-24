function [unew] = BoundaryCondition(boundaryType,u,boundaryValue) %,algorithm
%BOUNDARYCONDITION  Sets boundary condition for vector unew
unew=u;
switch strcmpi(boundaryType(1),'p') % check if first char is p or P
    case 0 % false, no periodic boundary
      % for i=0:m %where m=number of BC <- add later
      if strcmpi(boundaryType(1),'d') % Dirichlet boundary condition
        unew(end) = boundaryValue(end); 
        unew(1) = boundaryValue(1); 
      elseif strcmpi(boundaryType(1),'n') % Neumann boundary condition
        unew(end) =  u(end) - boundaryValue(end); % (unew(end-m+i)=unew(end-m)
        unew(1) = u(1) - boundaryValue(1); % (unew(m+1-i)=unew(m+2-i)
      else
        print('incorrect BC specification');
      end
    case 1 % true, periodic boundary
        unew(1)=u(end-1);
        unew(end)=unew(2);
end

end


% Author: Tyler James Pierce
% Tylerpierce644@gmail.com
% Created: Feb 07 2014

% Update History:



% Notes for me for future updates
%   if boundary is periodic
%   then u(end)=u(2) and u(1)=u(end-1)
%   elseif dirichlet
%   u(1) = g_D1(t)
%   u(end) = g_D2(t)
%   elseif neumann
%   unew(1) = u(1) - g_N1(t) and unew(end) = u(end) - g_N2(t)

% if alg=O[dx^m,dt^n]
% number of BC=m
% for i=1:imax
%     jp1=1+mod(i,imax);
%     jm1=imax-mod(imax+1-i,imax);
% end
