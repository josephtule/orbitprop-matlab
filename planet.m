classdef planet
    properties
        name
        radius
        mu
        j2 = 0
        rhos = 0;
        atmos_rot = 0;
    end

    methods
        function pl = planet(name,radius,mu,j2,rhos,atmos_rot)
            pl.name = name;
            pl.radius = radius;
            pl.mu = mu;
            pl.j2 = j2;
            pl.rhos = rhos;
            pl.atmos_rot = atmos_rot;
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

