function [q]=bc(q,area,imax,extrap,gamma,gamma1)

i1=imax-1;
i2=imax;

if strcmp(extrap,'True')==1
    q(i2,1)=q(i1,1);
    q(i2,2)=q(i1,2);
    q(i2,3)=q(i1,3);
end
u1 = q(i1,2)/q(i1,1);
p1 = gamma1*(q(i1,3) - .5*q(i1,2)*u1)/area(i1,1);
p2 = 1.931/gamma;
rho1 = q(i1,1)/area(i1,1);
s = (rho1^gamma)/p1;

rho2 = (s*p2)^(1/gamma);

c22 = gamma*p2/rho2;
c2 = sqrt(abs(c22));
c12 = gamma*p1/rho1;
c1 = sqrt(abs(c12));
r1 = u1+(2.0*c1/gamma1);
r2 = r1-(4.0*c2/gamma1);

u2 = .5*(r1+r2);
e2 = (p2/gamma1)+.5*rho2*u2*u2;

q(i2,1) = rho2*area(i2,1);
q(i2,2) = rho2*u2*area(i2,1);
q(i2,3) = e2*area(i2,1);


