function [u0] = StepWave(x,xStep)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

umin=0; %can use inputs for these later down the line
umax=1; %or an amplitude input
imax=numel(x);
u0=zeros(1,imax);
  for i=2:imax-1 % Generate initial wave  
      %if x(i)/x(end) < sqmin
    if x(i) < xStep 
      u0(i)=umin;
    else
      u0(i)=umax;
    end
  end

end

