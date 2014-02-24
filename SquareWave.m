function [u] = SquareWave(x,umin,umax)
%SQUAREWAVE Creates a square wave based upon input
%umax      ___   Square Wave
%umin____|    |______________

% configured parameters
sqmin=10; % 10 percent of total number of points
sqmax=20; % 20 percent of total number of points

imin=1;
imax=numel(x);

%imin=find(min(x)==x,1,'last'); %return the last index of x=xmin (ignores dummy points)
%imax=find(max(x)==x,1,'first'); %return the first index of x=xmax (ignores dummy points)

u=zeros(1,numel(x));
  for i=imin:imax % Generate initial wave  
    if i/imax < sqmin/100
      u(i)=umin;
    elseif i/imax <= sqmax/100
      u(i)=umax;
    else
      u(i)=umin;
    end
  end

end

