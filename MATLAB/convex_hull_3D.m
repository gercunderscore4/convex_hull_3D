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
		elseif abs(v(1)) >= abs(v(2)) && abs(v(1)) >= abs(v(3))
			if abs(v(2)) >= abs(v(3))
				u1 = [v(2) -v(1) 0];
			else
				u1 = [v(3) 0 -v(1)];
			end
		elseif abs(v(2)) >= abs(v(1)) && abs(v(2)) >= abs(v(3))
			if abs(v(1)) >= abs(v(3))
				u1 = [v(2) -v(1) 0];
			else
				u1 = [0 v(3) -v(2)];
			end
		else
			if abs(v(1)) >= abs(v(2))
				u1 = [v(3) 0 -v(1)];
			else
				u1 = [0 v(3) -v(2)];
			end
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
		v3_s = v3(indices,:);

		% inverse indices, 
		% so that they they can point to 
		% their counterparts on the unsorted lists
		[temp indexes] = sort(indices,'ascend');
		
		t_s = ones(len_p3,1);
		t_s(ii_s) = 0;
		
		disp(t_s(indexes)')
		
		% this is the tricky part
		
		% first, setup the counters
		% counting counter
		cc = 0;
		cc_max = len_p3+1;
		% j_previous
		jp = mod(cc-1-1,len_p3)+1;
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
		while cc <= cc_max && sum(t_s) > 2
			% get points
			pp = p3(indices(jp),:);
			pc = p3(indices(jc),:);
			pn = p3(indices(jn),:);
			
			t = dot( cross(pp-p,pn-p), pc-p );
			
			if t < 0
				t_s(jc) = 0;
				cc_max = cc + len_p3+1;
				jp = jp;
			else
				jp = jc;
			end
			jc = jn;
			% find next valid point
			cc = cc+1;
			jn = mod(cc-1,len_p3)+1;
			while t_s(jn) == 0
				cc = cc+1;
				jn = mod(cc-1,len_p3)+1;
			end
		end
		
		disp(t_s(indexes)')
		
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
		
		disp(t_s(indexes)')
		
	end
	
	% cull to proper size
	l3 = l3(1:len_l3,:);
end
