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
            'calc_coes',1, ...
            'animateopt',0, ...
            'frameskip',1, ...
            'solver','rk4', ...
            'plot3dopt',1, ...
            'plotcoesopt',1, ...
            'mission','orbit', ... % mission: orbit, launch, missile
            'date', [2018 1 17 16 20 36].' ...
            );

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
        animateopt
        frameskip
        solver
        plot3dopt
        plotcoesopt
        missilemode
        mission
        date
        specs = struct( ...
            'Cd',2.2, ...
            'area',(1e-3)^2/4, ...
            'mass', 100 ...
            );

    end

    methods
        function sc = spacecraft(config,specs)

            fn = fieldnames(config);
            for k = 1:numel(fn)
                sc.config.(fn{k}) = config.(fn{k});
            end
            pn = fieldnames(sc.config);
            for k = 1:numel(pn)
                sc.(pn{k}) = sc.config.(pn{k});
            end
            sn = fieldnames(specs);
            for k = 1:numel(sn)
                sc.specs.(sn{k}) = specs.(sn{k});
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

            if sc.plot3dopt == 1
                sc.plot3d();
            end
           
            if sc.plotcoesopt == 1
                sc.plotcoes();
            end

            if sc.animateopt == 1
                sc.drawvideo();
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
                    case "aero"
                        a_pert = a_pert + sc.aeropert(state);
                end
            end
            a = a + a_pert;
            state_dot(:,1) = [v;a];
        end

        function sc = proporbit(sc)

            N = ceil((sc.tspan(2) - sc.tspan(1))/sc.dt);
            
            for i = 1:N
                switch sc.solver
                    case 'rk4'
                        sc = rk4(sc,i);
                    case 'rk45'
                        sc = rk45(sc,i);
                end
                sc = proptime(sc,sc.date(:,i),i);
            end

        end

        function sc = rk4(sc,i)
            k1 = sc.dt * sc.eoms(sc.state(:,i));
            k2 = sc.dt * sc.eoms(sc.state(:,i)+k1/2);
            k3 = sc.dt * sc.eoms(sc.state(:,i)+k2/2);
            k4 = sc.dt * sc.eoms(sc.state(:,i)+k3);
            sc.state(:,i+1) = sc.state(:,i) + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
        end

        function sc = rk45(sc,i)
            %             A = [0 2/9 1/3 3/4 1 5/6];
            B = [0 0 0 0 0;
                2/9 0 0 0 0;
                1/12 1/4 0 0 0;
                69/128 -243/128 135/64 0 0;
                -17/12 27/4 -27/5 16/15 0;
                65/432 -5/16 13/16 4/27 5/144];
            %             C = [1/9 0 9/20 16/45 1/12 0];
            CH = [47/450 0 12/25 32/225 1/30 6/25];
            k1 = sc.dt * sc.eoms(sc.state(:,i));
            k2 = sc.dt * sc.eoms(sc.state(:,i) + B(2,1)*k1);
            k3 = sc.dt * sc.eoms(sc.state(:,i) + B(3,1)*k1 + B(3,2)*k2);
            k4 = sc.dt * sc.eoms(sc.state(:,i) + B(4,1)*k1 + B(4,2)*k2 + B(4,3)*k3);
            k5 = sc.dt * sc.eoms(sc.state(:,i) + B(5,1)*k1 + B(5,2)*k2 + B(5,3)*k3 + B(5,4)*k4);
            k6 = sc.dt * sc.eoms(sc.state(:,i) + B(6,1)*k1 + B(6,2)*k2 + B(6,3)*k3 + B(6,4)*k4 + B(6,5)*k5);
            sc.state(:,i+1) = sc.state(:,i) + CH(1)*k1 + CH(2)*k2 + CH(3)*k3 + CH(4)*k4 + CH(5)*k5 + CH(6)*k6;
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

        function sc = proptime(sc,date,i)
            sc.date(:,i+1) = date + sc.dt*[0;0;0;0;0;1];
            
            if sc.date(6,i+1) >= 60
                incnum = floor(sc.date(6,i+1)/60);
                sc.date(6,i+1) = mod(sc.date(6,i+1),60);
                sc.date(5,i+1) = sc.date(5,i+1) + incnum;
            end
            if sc.date(5,i+1) >= 60
                incnum = floor(sc.date(5,i+1)/60);
                sc.date(5,i+1) = mod(sc.date(5,i+1),60);
                sc.date(4,i+1) = sc.date(4,i+1) + incnum;
            end
            if sc.date(4,i+1) >= 24
                incnum = floor(sc.date(4,i+1)/24);
                sc.date(4,i+1) = mod(sc.date(4,i+1),24);
                sc.date(3,i+1) = sc.date(3,i+1) + incnum;
            end

            if ismember(sc.date(2,i+1),[2]) % if feb
                if mod(sc.date(1,i+1),4) == 0 % leapyear
                    maxdays = 29;
                else
                    maxdays = 28;
                end
            elseif ismember(sc.date(2,i+1),[1 3 5 7 8 10 12])
                maxdays = 31;
            elseif ismember(sc.date(2,i+1),[4 6 9 11])
                maxdays = 30;
            end
            
            if sc.date(3,i+1) > maxdays
                incnum = floor(sc.date(3,i+1)/maxdays);
                sc.date(3,i+1) = mod(sc.date(3,i+1),maxdays);
                sc.date(2,i+1) = sc.date(2,i+1) + incnum;
            end

            if sc.date(2,i+1) > 12
                incnum = floor(sc.date(2,i+1)/12);
                sc.date(2,i+1) = mod(sc.date(2,i+1),12);
                sc.date(1,i+1) = sc.date(1,i+1) + incnum;
            end
        end

        function a_pert = j2pert(sc,state)
            r = norm(state(1:3));
            a_pert(1,1) = -3/2 * sc.cb.j2 *  sc.cb.mu * sc.cb.radius^2 * state(1) / r^5 * (1 - 5*state(3)^2/r^2);
            a_pert(2,1) = -3/2 * sc.cb.j2 *  sc.cb.mu * sc.cb.radius^2 * state(2) / r^5 * (1 - 5*state(3)^2/r^2);
            a_pert(3,1) = -3/2 * sc.cb.j2 *  sc.cb.mu * sc.cb.radius^2 * state(3) / r^5 * (3 - 5*state(3)^2/r^2);
            % Vallado, 5ed, p597
        end

        function a_pert = aeropert(sc,state)
            alt = norm(state(1:3));
%             rho = atmos_density(alt);
            rho = 1;
            atm_rot = [1;1;1];
            v_rel = state(4:6) - cross(atm_rot,state(1:3));
            a_pert = -v_rel*.5*rho*norm(v_rel)*sc.specs.Cd*sc.specs.area/sc.specs.mass;
    
        end

        function atmos_density(sc,alt)
        end

        function find_rho(sc,alt)
        end

    

        function drawvideo(sc)
            xx = sc.state;
            set(0,'DefaultAxesFontName', 'Times New Roman')
            set(0,'DefaultAxesFontSize', 12)

            e = earth();
            line_width = 2;
            fontsize_labels = 14;


            figure(500)
            % Animate the motion
            set(gcf,'PaperPositionMode','auto')
            set(gcf, 'Color', 'w');
            set(gcf,'Units','normalized','OuterPosition',[0 0 0.55*0.7 1*0.7]);

            fileName = append('animations\animation-',datestr(datetime('now')),'.gif');
            fileName = strrep(fileName, ' ', '.');
            fileName = strrep(fileName, ':', '.');

            for k = 1:sc.frameskip:size(xx,2)
                e.plotplanet([0,0,0])

                hold on;
                x1 = xx(1,k); y1 = xx(2,k); z1 = xx(3,k);
                xtrack1 = xx(1,1:k); xtrack2 = xx(2,1:k); xtrack3 = xx(3,1:k);


                plot3(x1, y1, z1, 'k*', 'MarkerSize', 10); % plot position
                plot3(xtrack1,xtrack2,xtrack3);
                hold off
                ylabel('y-position','FontSize',fontsize_labels)
                xlabel('x-position','FontSize',fontsize_labels)
                zlabel('z-position','FontSize',fontsize_labels)

                % axis([0.985, 1.015 -0.010, 0.010]);
                axis equal
                box on;
                grid on

                ml = 10000;
                xlim([-ml, ml]);
                ylim([-ml, ml]);
                zlim([-ml, ml]);

                drawnow
                %lgd.FontSize = 7;
                %lgd.Location = 'eastoutside';
                exportgraphics(gcf,fileName,'Append',true); % save animation as GIF
                F(k) = getframe(gcf); % to get the current frame
            end
            close(gcf)
        end

        function plot3d(sc)
            figure()
            state0 = sc.state(:,1);
            plot3(sc.state(1,:),sc.state(2,:),sc.state(3,:))
            hold on
            plot3(state0(1),state0(2),state0(3),"*")
            quiver3(state0(1),state0(2),state0(3),state0(4)*1000,state0(5)*1000,state0(6)*1000)
            sc.cb.plotplanet([0,0,0])
            axis equal
            xlabel('x'),ylabel('y'),zlabel('z')
            grid on
        end

        function plotcoes(sc)
            figure(2)
            t = linspace(0,24*60*60 * 2,size(sc.state,2));
            titles = ["a","e","i","\Omega","\omega","\theta"];
            ylabels = ["km","e","degrees, \circ","degrees, \circ","degrees, \circ","degrees, \circ"];
            for k = 1:size(sc.coes,1)
                subplot(3,2,k)
                if ismember(k,[3,4,5,6])
                    plot(t,rad2deg(sc.coes(k,:)))
                else
                    plot(t,sc.coes(k,:))
                end
                ylabel(ylabels(k))
                title(titles(k))
                grid on
            end

        end
    end

end
