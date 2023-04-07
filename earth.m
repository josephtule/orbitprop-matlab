classdef earth < planet
    methods
        function e = earth()
            e@planet("earth", ... % name
                     6378.14, ... % radius
                     3.986004418e5, ... % mu
                     1.082635854e-3, ... % j2
                     [63.096 251.189 1000; 2.059e-4 5.909e-11 3.561e-15], ... % rhos, [altitude, density]
                     [0;0;72.9211e-6] ... % rotation vector
                     );
        end
    end
end
