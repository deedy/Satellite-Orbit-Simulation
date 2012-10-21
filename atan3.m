
function theta=atan3(x,y)
theta=atan2(y,x);
if (theta<0)
    theta=theta+2*pi;
end
