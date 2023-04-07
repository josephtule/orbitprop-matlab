classdef planet
    properties
        name
        radius
        mu
        j2
    end

    methods
        function pl = planet(name,radius,mu,j2)
            pl.name = name;
            pl.radius = radius;
            pl.mu = mu;
            pl.j2 = j2;
        end

        function plotplanet(pl,offset)
            [x,y,z] = sphere;
            x = x*pl.radius+offset(1);
            y = y*pl.radius+offset(2);
            z = z*pl.radius+offset(3);
            surf(x,y,z)
        end
    end
end

