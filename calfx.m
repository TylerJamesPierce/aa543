function [flux]=calfx(q,flux,istart,imax,gamma1)
 for i=istart:imax 
    u=q(i,2)/q(i,1);
    pa=gamma1*(q(i,3)-.5*q(i,2)*u);
    flux(i,1)=q(i,2);
    flux(i,2)=q(i,2)*u + pa;
    flux(i,3)=u*(q(i,3)+pa);
 end
