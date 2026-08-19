 function [Vmag, Vx, Vy, Gx, Gy, meta] = circle_method_cv(T, dx, dy, r)
            % Circle Method (Siles-Paredes et al., 2022) local CV from a LAT map.
            %   T     : activation-time map (ms); NaN/0 = off-tissue
            %   dx,dy : pixel size (mm);  r : circle radius (pixels)
            % Returns CV magnitude in mm/ms (= m/s) and direction-bearing components
            % Vx,Vy — same convention/units as conduction_velocity.m, so the caller's
            % "*100" yields cm/s and atan2d(Vy,Vx) gives propagation direction.
            try
                if nargin < 4 || isempty(r), r = 10; end
                r = max(4, round(r));
                dtheta = 60;  nang = 36;  speed_max = 2.50;   % mm/ms cap (= 250 cm/s)
                ps = (dx + dy)/2;                              % mm per pixel
                [nr, nc] = size(T);
                T = double(T);  T(T == 0) = NaN;
                validT = isfinite(T);
                Tfill  = T;  Tfill(~validT) = 0;
                mask   = double(validT);
                [cols, rows] = meshgrid(1:nc, 1:nr);          % X=col, Y=row
                angles = linspace(0, 180, nang+1);  angles(end) = [];

                dLAT = nan(nr, nc, nang);               % SIGNED ΔLAT per orientation
                for k = 1:nang
                    th  = angles(k);
                    dxr = r*cosd(th);  dyr = r*sind(th);      % chord half-vector (col,row)
                    lp = interp2(cols, rows, Tfill, cols + dxr, rows + dyr, 'linear', NaN);
                    lm = interp2(cols, rows, Tfill, cols - dxr, rows - dyr, 'linear', NaN);
                    mp = interp2(cols, rows, mask,  cols + dxr, rows + dyr, 'linear', 0);
                    mm = interp2(cols, rows, mask,  cols - dxr, rows - dyr, 'linear', 0);
                    d  = lp - lm;                             % signed: + => later on the +th side
                    d(mp < 0.999 | mm < 0.999) = NaN;        % chord left the tissue
                    dLAT(:,:,k) = d;
                end
                absd = abs(dLAT);

                % propagation direction = chord with the largest |ΔLAT|
                [maxd, pdidx] = max(absd, [], 3, 'omitnan');
                th_pd = angles(pdidx);

                % resolve the 180-deg ambiguity: orient toward LATER activation
                % (early->late = away from origin), matching the gradient convention
                [iir, jjc] = ndgrid(1:nr, 1:nc);
                good_pd = isfinite(maxd);
                lin = sub2ind([nr, nc, nang], iir(good_pd), jjc(good_pd), pdidx(good_pd));
                s   = sign(dLAT(lin));                        % +1 if +th_pd is the later side
                thv = th_pd(good_pd);
                thv(s < 0) = thv(s < 0) + 180;
                th_pd(good_pd) = thv;
                th_pd(~good_pd) = NaN;

                % cosine-weighted, back-projected speed over the wedge around th_pd
                acc = zeros(nr, nc);  cnt = zeros(nr, nc);
                for k = 1:nang
                    delta  = abs(mod(angles(k) - th_pd + 90, 180) - 90);   % dist to PD
                    d      = absd(:,:,k);
                    spd    = (2*r*ps) .* cosd(delta) ./ d;                 % mm/ms
                    good   = (delta <= dtheta/2) & isfinite(spd) & (d > 1e-6);
                    acc(good) = acc(good) + spd(good);
                    cnt(good) = cnt(good) + 1;
                end
                Vmag = nan(nr, nc);
                Vmag(cnt > 0) = acc(cnt > 0) ./ cnt(cnt > 0);
                Vmag(Vmag > speed_max) = NaN;   % reject flat-LAT blow-ups (don't clamp to the cap)
                Vmag(~validT) = NaN;

                Vx = Vmag .* cosd(th_pd);     % direction-bearing (axis in [0,180))
                Vy = Vmag .* sind(th_pd);
                Gx = cosd(th_pd) ./ Vmag;     % |∇T| components (ms/mm), for display
                Gy = sind(th_pd) ./ Vmag;

                meta = struct('method','circle','radius',r,'dtheta',dtheta, ...
                    'nang',nang,'dx',dx,'dy',dy,'validT',validT, ...
                    'good',isfinite(Vmag));
            catch ME
                Vmag = nan(size(T)); Vx = Vmag; Vy = Vmag; Gx = Vmag; Gy = Vmag;
                meta = struct('method','circle','error',ME.message,'good',false(size(T)));
                
            end
        end
