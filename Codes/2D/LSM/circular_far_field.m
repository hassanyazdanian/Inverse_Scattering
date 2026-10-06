function F = circular_far_field(k,n,radius,angles,sound_soft,center)
%CIRCULAR_FAR_FIELD Analytical disk far field, observation rows/incident columns.
% Convention: us = exp(ikr)/sqrt(r) * u_infinity + O(r^(-3/2)).
if nargin < 5, sound_soft = false; end
if nargin < 6, center = [0;0]; end
validateattributes(k,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(n,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(radius,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(angles,{'numeric'},{'vector','real','finite','nonempty'});
validateattributes(center,{'numeric'},{'vector','numel',2,'real','finite'});
angles = angles(:); center = center(:);
delta = bsxfun(@minus,angles,angles');
F = complex(zeros(numel(angles)));
x = k*radius; y = sqrt(n)*x;
max_order = ceil(max(x,y)+15);
for m = -max_order:max_order
    J = besselj(m,x); H = besselh(m,1,x);
    if sound_soft
        coefficient = -J/H;
    else
        Jp = (besselj(m-1,x)-besselj(m+1,x))/2;
        Hp = (besselh(m-1,1,x)-besselh(m+1,1,x))/2;
        Ji = besselj(m,y); Jip = sqrt(n)*(besselj(m-1,y)-besselj(m+1,y))/2;
        coefficient = (J*Jip-Jp*Ji)/(Ji*Hp-Jip*H);
    end
    F = F+coefficient*exp(1i*m*delta);
end
projection = [cos(angles) sin(angles)]*center;
F = sqrt(2/(pi*k))*exp(-1i*pi/4)*F.*exp(1i*k*bsxfun(@minus,projection',projection));
end
