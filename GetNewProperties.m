function [unew] = GetNewProperties(cfl,u0) %algorithm, iter, du, unew, uold, boundary, theta, imax
%GETNEWPROPERTIES  Updates the values of the property of u, based on the
%Courant Number, cfl, and the initial properties u, (and eventually alg, current iteration, etc)

du=zeros(size(u0));

% if algorithm==Lax-Wendroff 
  for i = 2:(numel(u0)-1) % update all values except boundary values
%     up = ( u0(i+1) + u0(i) ) / 2 - cfl * (u0(i+1) - u0(i)); % calc property at u(i+1/2)
%     um = ( u0(i) + u0(i-1) ) / 2 - cfl * (u0(i) - u0(i-1)); % calc property at u(i-1/2)
%     du(i) = - cfl * ( up - um ); % combine with corrector step
    du(i)= cfl^2/2*(u0(i+1)+u0(i-1)-2*u0(i))-(cfl/2)*(u0(i+1)-u0(i-1));
  end
  unew = u0 + du; % update u

end

% Author: Tyler James Pierce
% Tylerpierce644@gmail.com
% Created: Feb 07 2014

% Update History:
% 02/17/2014 - Changed Lax-Wendroff to a two-step method
