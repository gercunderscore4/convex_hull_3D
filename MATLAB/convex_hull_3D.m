%{
	3D convex hull
	assumes p3 is a convex set
	returns list of triangles, size = n x 3
%}
function l3 = convex_hull_3D(p3)
	len_p3 = size(p3,1);

	% maximum possible number of triangles (using this algorithm)
	l3 = zeros(len_p3^2,3);
	len_l3 = 0;

	% find center
	c3 = [min(p3(:,1))+max(p3(:,1)) min(p3(:,2))+max(p3(:,2)) min(p3(:,3))+max(p3(:,3))]/2;
	
	% outward vectors
	v3 = bsxfun(@minus, p3, c3);
	
	ii = 0;
	ii = find(p3(:,1)==max(p3(:,1)),1);

	% cs marks which points form the convex set
	% 0 means not registered as convex (yet)
	% 1 means registered, but not acted upon
	% 2 means registered and acted upon
	cs = zeros(size(p3,1),1);
	cs(ii) = 1;

	while ii ~= 0
		
		% for ease and fewer look-ups
		p = p3(ii,:);
		v = v3(ii,:);
		
		% get directional vectors
		% don't bother with the case v3(ii,:) == [0 0 0]
		% if that happens, something already went wrong
		if v(1) == 0
			u1 = [1 0 0];
		elseif v(2) == 0
			u1 = [0 1 0];
		else
			u1 = [v(2) -v(1) 0];
		end
		u2 = cross(v,u1);
		
		% get angles of points around v
		theta = zeros(len_p3,1);
		for jj = 1:len_p3
			theta(jj) = atan2( dot(u2,v3(jj,:)), dot(u1,v3(jj,:)) );
		end
		
		% sort ccw
		[theta indices] = sort(theta, 'ascend');
		ii_s = find(indices == ii, 1);
		% fan test list
		t_s = ones(len_p3,1);
		t_s(ii_s) = 0;
		
		% DEBUG
		[temp indexes] = sort(indices,'ascend');
		
		% DEBUG
		%disp(t_s(indexes)')
		
		% this is the tricky part
		
		% first, setup the counters
		% j0
		% j1
		% j2
		% j3 == cc
		% counting counter
		cc = 0;
		cc_max = len_p3;
		% j_back
		j1 = imod(cc-1,len_p3);
		while t_s(j1) == 0
			j1 = imod(j1-1,len_p3);
		end
		% j_previous
		j0 = imod(j1-1,len_p3);
		while t_s(j0) == 0
			j0 = imod(j0-1,len_p3);
		end
		% j_current
		j2 = imod(cc,len_p3);
		while t_s(j2) == 0
			cc = cc+1;
			j2 = imod(cc,len_p3);
		end
		% j_next
		% cc tracks this one because it's the leading edge
		cc = cc+1;
		j3 = imod(cc,len_p3);
		while t_s(j3) == 0
			cc = cc+1;
			j3 = imod(cc,len_p3);
		end
		while cc <= cc_max
			
			if dot( cross(p3(indices(j1),:)-p,p3(indices(j3),:)-p), p3(indices(j2),:)-p ) < 0
				t_s(j2) = 0;
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				
				% move ahead
				j2 = j3;
				% find next valid point
				cc = cc+1;
				j3 = imod(cc,len_p3);
				while t_s(j3) == 0
					cc = cc+1;
					j3 = imod(cc,len_p3);
				end
				
				t_b = 0;
				t_c = 0;
				while t_b == 0 || t_c == 0
					
					% check backwards
					t_s(j1) = 1*(  dot( cross(p3(indices(j0),:)-p,p3(indices(j2),:)-p), p3(indices(j1),:)-p ) >= 0  );
					t_b = t_s(j1);
					% if failed
					if t_b == 0
						% move back
						j1 = j0;
						j0 = imod(j0-1,len_p3);
						while t_s(j0) == 0
							j0 = imod(j0-1,len_p3);
						end
					end

					% check forwards
					t_s(j2) = 1*(  dot( cross(p3(indices(j1),:)-p,p3(indices(j3),:)-p), p3(indices(j2),:)-p ) >= 0  );
					t_c = t_s(j2);
					% if failed
					if t_c == 0
						% move ahead
						j2 = j3;
						% find next valid point
						cc = cc+1;
						j3 = imod(cc,len_p3);
						while t_s(j3) == 0
							cc = cc+1;
							j3 = imod(cc,len_p3);
						end
					end

				end
				
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			else
				%t_s(j2) = 1;
				
				j0 = j1;
				j1 = j2;
				j2 = j3;
				% find next valid point
				cc = cc+1;
				j3 = imod(cc,len_p3);
				while t_s(j3) == 0
					cc = cc+1;
					j3 = imod(cc,len_p3);
				end
			end
			
		end
		
		% DEBUG
		%disp(t_s(indexes)')
				
		% create triangles
		
		% counters
		cc = 0;
		cc_max = len_p3-1;
		% j_previous
		j0 = imod(cc-1,len_p3);
		while t_s(j0) == 0
			j0 = imod(j0-1,len_p3);
		end
		while cc <= cc_max
			
			j2 = imod(cc,len_p3);
			if t_s(j2) == 1
				% record triangle
				i2 = indices(j2);
				i3 = indices(j0);
			
				if i2 < ii || i3 < ii
					% already accounted for
					%t_s(j2) = 0;
				elseif i2 < i3
					len_l3 = len_l3+1;
					l3(len_l3,:) = [ii i2 i3];
				else % if i2 < i3
					len_l3 = len_l3+1;
					l3(len_l3,:) = [ii i3 i2];
				end
				%
				
				j0 = j2;
			end
			cc = cc+1;
		end
		
		% DEBUG
		%disp(t_s(indexes)')
		
		cs(ii) = 2;
		ii = 0;
		for jj = 1:len_p3
			if t_s(indexes(jj)) == 1 && cs(jj) == 0
				cs(jj) = 1;
			end
			if cs(jj) == 1
				ii = jj;
			end
		end

	end
	
	% cull to proper size
	l3 = l3(1:len_l3,:);
end
