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
            'mass',[], ...
            'calc_coes',0);

        cb
        frame
        dt
        tspan
        coes
        state
        perts
        stopcon
        mass
        calc_coes

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

            sc = sc.centralbody();
            if isempty(sc.state)
                sc = sc.coes2state();
            end

            sc = sc.proporbit();

            if sc.calc_coes == 1
                for i = 1:size(sc.state,2)-1
                    sc.coes(:,i+1) = sc.state2coes(sc.state(:,i+1));
                end
            end


        end

        function state_dot = eoms(sc,state)
            r(1:3,1) = state(1:3);
            v(1:3,1) = state(4:6);
            a(1:3,1) = -r * sc.cb.mu / norm(r)^3;
            a_pert = zeros(3,1);
            for i = 1:length(sc.perts)
                switch sc.perts(i)
                    case "j2"
                        a_pert = a_pert + sc.j2pert(state);
                    case "srp"
                    case "thirdbody"
                end
            end
            a = a + a_pert;
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
            a = sc.coes(1); ecc = sc.coes(2); i = sc.coes(3);
            RAAN = sc.coes(4); AOP = sc.coes(5); TA = sc.coes(6);

            p = a*(1-ecc^2);
            xyz = [p*cos(TA)/(1+ecc*cos(TA)); p*sin(TA)/(1+ecc*cos(TA)); 0];
            xyz_dot = [-sqrt(sc.cb.mu/p)*sin(TA); sqrt(sc.cb.mu/p)*(ecc+cos(TA)); 0];

            NP = sc.rot_N2P(); % from inertial to perifocal
            PN = NP'; % from perifocal to inertial
            sc.state(1:3,1) = PN*xyz;
            sc.state(4:6,1) = PN*xyz_dot;
        end

        function cc = state2coes(sc,state)
            r_vec = state(1:3); v_vec = state(4:6); h_vec = cross(r_vec,v_vec);
            r = norm(r_vec); v = norm(v_vec); h = norm(h_vec);
            energy = v^2/2 - sc.cb.mu/r;
            a = -sc.cb.mu/2/energy;
            e_vec = cross(v_vec,h_vec)/sc.cb.mu - r_vec/r;
            ecc = norm(e_vec);
            k_hat = [0;0;1]; n_vec = cross(k_hat,h_vec);
            
            i_e = e_vec/ecc; i_h = h_vec/h; i_p = cross(i_h,i_e);

            R_NP = [i_e.';i_p.';i_h.']; % from N to P
            h_inert = R_NP' * h_vec;
            n_inert = R_NP' * n_vec;
            e_inert = R_NP' * e_vec;

            i = acos(R_NP(3,3));
            RAAN = atan2(R_NP(3,1),-R_NP(3,2));
            AOP = atan2(R_NP(1,3),R_NP(2,3));
            TA = acos(dot(e_vec,r_vec)/norm(e_vec)/norm(r_vec));
            if n_inert(2) < 0
                RAAN = 2*pi-RAAN;
            end
            if e_inert(3) < 0
                AOP  = 2*pi - AOP;
            end
            if dot(r_vec,v_vec) < 0
                TA = 2*pi - TA;
            end
            cc = [a,ecc,i,RAAN,AOP,TA].';
        end


        function R = rot_N2P(sc)
            % rotation matrix from the GCI (inertial) frame to Perifocal
            % frame
            i = sc.coes(3); RAAN = sc.coes(4); AOP = sc.coes(5);
            R(1,1) = cos(RAAN)*cos(AOP)-sin(RAAN)*sin(AOP)*cos(i);
            R(1,2) = sin(RAAN)*cos(AOP)+cos(RAAN)*sin(AOP)*cos(i);
            R(1,3) = sin(AOP)*sin(i);
            R(2,1) = -cos(RAAN)*sin(AOP)-sin(RAAN)*cos(AOP)*cos(i);
            R(2,2) = -sin(RAAN)*sin(AOP)+cos(RAAN)*cos(AOP)*cos(i);
            R(2,3) = cos(AOP)*sin(i);
            R(3,1) = sin(RAAN)*sin(i);
            R(3,2) = -cos(RAAN)*sin(i);
            R(3,3) = cos(i);
        end

        function a_pert = j2pert(sc,state)
            r = norm(state(1:3));
            a_pert(1,1) = -3/2 * sc.cb.j2 *  sc.cb.mu * sc.cb.radius^2 * state(1) / r^5 * (1 - 5*state(3)^2/r^2);
            a_pert(2,1) = -3/2 * sc.cb.j2 *  sc.cb.mu * sc.cb.radius^2 * state(2) / r^5 * (1 - 5*state(3)^2/r^2);
            a_pert(3,1) = -3/2 * sc.cb.j2 *  sc.cb.mu * sc.cb.radius^2 * state(3) / r^5 * (3 - 5*state(3)^2/r^2);
            % Vallado, 5ed, p597
        end
    end

end
