function blast(x,y,z,R)
    for i = 0:0.5:2*pi
        for j=0:0.5:2*pi
                plot3 (x+R*cos(j)*cos(i), y+R*cos(j)*sin(i), z+R*sin(j), '-r*');
        end
    end
end