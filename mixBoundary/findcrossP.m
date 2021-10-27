function P = findcrossP(K,vec)
vec = vec/norm(vec);
sx1 = (1-K(1))/vec(1);
sx2 = (-K(1))/vec(1);
if abs(sx1)+abs(sx2)<= sqrt(3)
    if abs(sx1)<abs(sx2)
        sx = sx1;
    else
        sx = sx2;
        
    end
else
    sx = 0;
end

sy1 = (1-K(2))/vec(2);
sy2 = (-K(2))/vec(2);
if abs(sy1)+abs(sy2)<=sqrt(3)
    if abs(sy1)<abs(sy2)
        sy = sy1;
    else
        sy = sy2;
        
    end
else
    sy = 0;
end

sz1 = (1-K(3))/vec(3);
sz2 = (-K(3))/vec(3);
if abs(sz1)+abs(sz2)<=sqrt(3)
    if abs(sz1)<abs(sz2)
        sz = sz1;
    else
        sz = sz2;
        
    end
else
    sz = 0;
end
sx
sy
sz
if sx ~= 0
    s = sx;
elseif sy ~= 0
    s = sy;
else
    s = sz;
end
P = K + s*vec;
