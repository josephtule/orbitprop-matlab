classdef earth < planet
    methods
        function e = earth()
            e@planet("earth", ... % name
                     6378.14, ... % radius
                     3.986004418e5, ... % mu
                     1.082635854e-3); % j2
        end
    end
end
