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
	
	for ii = 1:len_p3
		
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
		%[temp indexes] = sort(indices,'ascend');
		
		% DEBUG
		%disp(t_s(indexes)')
		
		% this is the tricky part
		
		% first, setup the counters
		% jp
		% jb
		% jc
		% jn == cc
		% counting counter
		cc = 0;
		cc_max = len_p3;
		% j_back
		jb = mod(cc-1-1,len_p3)+1;
		while t_s(jb) == 0
			jb = mod(jb-1-1,len_p3)+1;
		end
		% j_previous
		jp = mod(jb-1-1,len_p3)+1;
		while t_s(jp) == 0
			jp = mod(jp-1-1,len_p3)+1;
		end
		% j_current
		jc = mod(cc-1,len_p3)+1;
		while t_s(jj) == 0
			cc = cc+1;
			jj = mod(cc-1,len_p3)+1;
		end
		% j_next
		% cc tracks this one because it's the leading edge
		cc = cc+1;
		jn = mod(cc-1,len_p3)+1;
		while t_s(jn) == 0
			cc = cc+1;
			jn = mod(cc-1,len_p3)+1;
		end
		while cc <= cc_max
			
			if dot( cross(p3(indices(jb),:)-p,p3(indices(jn),:)-p), p3(indices(jc),:)-p ) < 0
				t_s(jc) = 0;
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				
				% move ahead
				jc = jn;
				% find next valid point
				cc = cc+1;
				jn = mod(cc-1,len_p3)+1;
				while t_s(jn) == 0
					cc = cc+1;
					jn = mod(cc-1,len_p3)+1;
				end
				
				t_b = 0;
				t_c = 0;
				while t_b == 0 || t_c == 0
					
					% check backwards
					t_s(jb) = 1*(  dot( cross(p3(indices(jp),:)-p,p3(indices(jc),:)-p), p3(indices(jb),:)-p ) >= 0  );
					t_b = t_s(jb);
					% if failed
					if t_b == 0
						% move back
						jb = jp;
						jp = mod(jp-1-1,len_p3)+1;
						while t_s(jp) == 0
							jp = mod(jp-1-1,len_p3)+1;
						end
					end

					% check forwards
					t_s(jc) = 1*(  dot( cross(p3(indices(jb),:)-p,p3(indices(jn),:)-p), p3(indices(jc),:)-p ) >= 0  );
					t_c = t_s(jc);
					% if failed
					if t_c == 0
						% move ahead
						jc = jn;
						% find next valid point
						cc = cc+1;
						jn = mod(cc-1,len_p3)+1;
						while t_s(jn) == 0
							cc = cc+1;
							jn = mod(cc-1,len_p3)+1;
						end
					end

				end
				
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			else
				%t_s(jc) = 1;
				
				jp = jb;
				jb = jc;
				jc = jn;
				% find next valid point
				cc = cc+1;
				jn = mod(cc-1,len_p3)+1;
				while t_s(jn) == 0
					cc = cc+1;
					jn = mod(cc-1,len_p3)+1;
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
		jp = mod(cc-1-1,len_p3)+1;
		while t_s(jp) == 0
			jp = mod(jp-1-1,len_p3)+1;
		end
		while cc <= cc_max
			
			jc = mod(cc-1,len_p3)+1;
			if t_s(jc) == 1
				% record triangle
				i2 = indices(jc);
				i3 = indices(jp);
			
				if i2 < ii || i3 < ii
					% already accounted for
					t_s(jc) = 0;
				elseif i2 < i3
					len_l3 = len_l3+1;
					l3(len_l3,:) = [ii i2 i3];
				else % if i2 < i3
					len_l3 = len_l3+1;
					l3(len_l3,:) = [ii i3 i2];
				end
				%
				
				jp = jc;
			end
			cc = cc+1;
		end
		
		% DEBUG
		%disp(t_s(indexes)')
		
	end
	
	% cull to proper size
	l3 = l3(1:len_l3,:);
end
