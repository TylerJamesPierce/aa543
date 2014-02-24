function [unew] = GaussianWave(range,x,xpeak) %dx
%GAUSSIANWAVE returns a gauss wave initial condition noramlized
xmin=min(range);
xmax=max(range);
imax=numel(x);
unew=zeros(1,imax);
xwidth=x(3)-x(2);
% or
%xt(i)=mod(x(i)-t,pi);
% Gaussian Wave
%   uex(i)=1+exp(-50*(xt(i)-0.5)^2);
%can use dx(i) or x(i)-x(i-1)

for i=1:imax
    unew(i) = 1+exp(-2/(xwidth)*((x(i)-xpeak)*(xmax+xmin)/(xmax-xmin))^2);
end

end


