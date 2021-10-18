classdef tomo 
    methods (Static)
        function f = diskChar(a,b,h,k,r)
            dist = sqrt((a-h)^2+(b-k)^2);
            if dist <= r
                f = 1;
            end 
            if dist > r
                f = 0;
            end 
        end 
        function F = F(phi, t, d)
            M_phi = [cos(phi) -sin(phi); sin(phi), cos(phi)];
            a = d/2;
            b = sqrt(d^2-1)/2;
            gamma_t = [a*cos(t)+.5; b*sin(t)];
            g_p = mtimes(M_phi, gamma_t);
            centerA = 0;
            centerB = 0;
            radius = 1/2;
            F = tomo.diskChar(g_p(1), g_p(2), centerA, centerB, radius)*sqrt(d^2-cos(t)^2)/2;
        end 
        function trap = trap(phi, d, l)
            delta_t = (2*pi)/(2*l);
            sum = 0;
            for i = 1:2*l-1
                sum = sum + tomo.F(phi, delta_t*i, d);
            end
            trap = (delta_t/2)*(tomo.F(phi, 0, d)+tomo.F(phi, 2*pi, d)+2*sum);
        end
        function Rf = Rf(phi_a, phi_b, k,l)
            % k = intervals for phi 
            % l = intervals for d and for t 
            % IS T determined by delta_phi 
            d_max = 7;
            delta_d = (d_max-1)/l;
            delta_phi = (phi_a-phi_b)/k;
            Rf = [];
            for i = 0:l-1
                d = 1 + delta_d*i;
                for j = 0:k-1
                    phi = delta_phi*j;
                    A_element = tomo.trap(phi, d, l);
                    Rf(i+1, j+1) = A_element;
                end
            end
            Rf_smooth = [];
            for i = 1:l
                for j = 1:k
                    if i == 1
                        Rf_smooth(i,j) = 1/6*Rf(3,j)+2/6*Rf(2,j)+3/6*Rf(1,j);
                    end
                    if i == 2
                        Rf_smooth(i,j) = 1/8*Rf(4,j)+2/8*Rf(3,j)+3/8*Rf(2,j)+2/8*Rf(1,j);
                    end
                    if i == l
                        Rf_smooth(i,j) = 1/6*Rf(l-2,j)+2/6*Rf(l-1,j)+3/6*Rf(l,j);
                    end
                    if i == l-1
                        Rf_smooth(i,j) = 1/8*Rf(l-3,j)+2/8*Rf(l-2,j)+3/8*Rf(l-1,j)+2/8*Rf(l,j);
                    end
                    if i ~= 1 && i ~= 2 && i ~= l && i ~= l-1 
                        Rf_smooth(i,j) = 1/9*Rf(i-2,j)+2/9*Rf(i-1,j)+3/9*Rf(i,j)+2/9*Rf(i+1,j)+1/9*Rf(i+2,j);
                    end 
                end 
            end            
        end
        function der = derived(delta_d, Rf,l,k)
            der = [];
            for i = 1:l
                if i == 1
                    for j = 0:k-1
                        der(1, j+1) = (Rf(1,j+1)-2*Rf(2, j+1) + Rf(3, j+1))/(delta_d^2);
                    end
                end 
                if i == l
                    for j = 0:k-1
                        der(i, j+1) = (Rf(i-2,j+1)-2*Rf(i-1, j+1)+ Rf(i, j+1))/(delta_d^2);
                    end 
                end
                if i ~= 1 && i~= l
                    for j = 0:k-1
                        der(i,j+1) = (Rf(i-1, j+1)-2*Rf(i,j+1)+Rf(i+1,j+1))/(delta_d^2);
                    end
                end 
            end 
        end 
        
        
        function d = d(phi, x1, x2)
            norm = sqrt(x1^2+x2^2);
            d = norm + sqrt((x1-cos(phi))^2+(x2-sin(phi))^2);
        end 
        
        function interpolation = interpolation(j,p, n, data, d)
            %Change D boundary currently [1,5]
            delta_d = 4/n;
            dj1 = 1+ (j-1)*delta_d;
            j2 = j-1;
            interpolation = (data(j,p)-data(j2, p))*((d-dj1)/delta_d)+ data(j2,p);
        end
        
        function BP = BP(n,k,l)
            %change dmax 
            dmax = 7;
            delta_d = (dmax-1)/l;
            %change x_max 
            xmax = 2;
            delta_x = (2*xmax)/n;
            %change phi bounds 
            phi_a = 0; 
            phi_b = 2*pi;
            delta_phi = (phi_a-phi_b)/k;
            Rf = tomo.Rf(phi_a, phi_b, k,l);
            derived = tomo.derived(delta_d,Rf,l,k); 
            BP_Rf =[];
            BP_derived = [];
            for i= 0:n-1
                x2 = xmax-i*(delta_x);
                for j = 0:n-1
                    x1 = -xmax+j*delta_x;
                    Rf_trap = 0;
                    derived_trap = 0;
                    for p = 0:k-1 
                        phi = phi_a+delta_phi*p;
                        d = tomo.d(phi, x1, x2);
                        j_val = floor(((d-1)*n)/(dmax-1))+2;
                        if j_val == 1
                            j_val = 2;
                        end 
                        Rf_interpolation = tomo.interpolation(j_val, p+1, n, Rf, d);
                        derived_interpolation = tomo.interpolation(j_val, p+1, n, derived, d);
                        if p == 0 || p == k-1
                            Rf_trap = Rf_trap + Rf_interpolation;
                            derived_trap = derived_trap + derived_interpolation;
                        end 
                        if p ~=0 && p~= k-1
                            Rf_trap = Rf_trap + 2*Rf_interpolation; 
                            derived_trap = derived_trap + 2*derived_interpolation;
                        end 
                        BP_Rf(i+1, j+1) = (delta_phi/2)*Rf_trap; 
                        BP_derived(i+1, j+1) = (delta_phi/2)*derived_trap;
                    end 
                end 
            end
            centerA = 0;
            centerB = 0;
            radius = 1/2;
            shapeRadius = radius*l/4;
            shapeCenterA = centerA*(l/4)+l/2; 
            shapeCenterB = l-(centerB*(l/4)+l/2);
        %   figure(1); imagesc(BP_Rf);
            figure(1); imagesc(BP_Rf);
            hold on 
         %   p0 = nsidedpoly(1000, 'Center', [shapeCenterA, shapeCenterB], 'Radius', shapeRadius);
          %  plot(p0, 'FaceColor', 'r', 'FaceAlpha',.01)
           % p00 = nsidedpoly(1000, 'Center', [l/2, l/2], 'Radius', l/4);
           % plot(p00, 'FaceColor', 'b', 'FaceAlpha',.01)
            p000 = nsidedpoly(1000, 'Center', [shapeCenterA, shapeCenterB], 'Radius', shapeRadius);
        %    plot(p000, 'FaceColor', 'r', 'FaceAlpha',0)
            p0000 = nsidedpoly(1000, 'Center', [l/2, l/2], 'Radius', l/4);
            xlabel('x_1 values')
            ylabel('x_2 values')
          %  plot(p0000, 'FaceColor', 'b', 'FaceAlpha',0)
            hold off 
            axis square; 
          %  figure(2); plot(Rf(:, 200)); 
            figure(2); imagesc(BP_derived);
            hold on 
            p0 = nsidedpoly(1000, 'Center', [shapeCenterA, shapeCenterB], 'Radius', shapeRadius);
           % plot(p0, 'FaceColor', 'r', 'FaceAlpha',.1, 'LineStyle','none')
            p00 = nsidedpoly(1000, 'Center', [l/2, l/2], 'Radius', l/4);
          %  plot(p00, 'FaceColor', 'b', 'FaceAlpha',.1, 'LineStyle','none')
            xlabel('x_1 values')
            ylabel('x_2 values')
            hold off 
            axis square;
            BP = 0;
        end 
    end
end
            