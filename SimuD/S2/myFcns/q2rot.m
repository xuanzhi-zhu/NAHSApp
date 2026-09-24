function rot = q2rot(q)
S=[0,-1;1,0];

rot=eye(2) + 2.*q(1).*q(2).*S + 2.*q(2).*q(2).*S*S;

end

