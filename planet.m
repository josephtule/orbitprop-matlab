classdef spacecraft
    properties
        % default config values
        config = struct( ...
            'cb', 'earth', ...
            'frame', 'J2000', ...
            'dt', [], ...
            'tspan', [],...
            'coes', [], ...
            'state', [], ...
            'perts', [], ...
            'stopcon',[], ...
            'mass',[]);
       
        cb
        frame
        dt
        tspan
        coes
        state
        perts
        stopcon
        mass

    end


    methods 
        function sc = spacecraft(config)

            fn = fieldnames(config);
            for k = 1:numel(fn)
                sc.config.(fn{k}) = config.(fn{k});
            end
            pn = fieldnames(sc.config);
            for k = 1:numel(pn)
                sc.(pn{k}) = sc.config.(pn{k});
            end
            
            % if isempty(state)
            %     sc = sc.coes2state();
            % end



            sc = sc.centralbody();
            sc = sc.proporbit();
        end

        function state_dot = eoms(sc,state)
            r(1:3,1) = state(1:3);
            v(1:3,1) = state(4:6);
            a(1:3,1) = -r * sc.cb.mu / norm(r)^3;
%             for i = 1:length(perts)
%                 a = a + a_pert;
%             end
            state_dot(:,1) = [v;a];
        end

        function sc = proporbit(sc)

            N = ceil((sc.tspan(2) - sc.tspan(1))/sc.dt);
    
            for i = 1:N
                k1 = sc.dt * sc.eoms(sc.state(:,i));
                k2 = sc.dt * sc.eoms(sc.state(:,i)+k1/2);
                k3 = sc.dt * sc.eoms(sc.state(:,i)+k2/2);
                k4 = sc.dt * sc.eoms(sc.state(:,i)+k3);
                sc.state(:,i+1) = sc.state(:,i) + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
            end
        end
        
        function sc = centralbody(sc)
            switch sc.cb
                case 'earth'
                    sc.cb = earth();
                otherwise
                    sc.cb = earth();
            end
        end

        function sc = coes2state(sc)
            
        
        end

    end

end
