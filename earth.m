classdef earth < planet
    methods
        function e = earth()
            e@planet("earth",6378.14,3.986004418e5);
        end
    end
end
