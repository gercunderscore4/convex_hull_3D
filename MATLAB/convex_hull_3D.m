%{
	3D convex hull
	assumes p3 is a convex set
	returns list of triangles, size = n x 3
	
	Create a list of vectors from the center of the object.
	For every point/vector (pv):
		Create the unit vectors, u1 and u2, perpendicular to u3.
		Use u1 and u2 to sort the points ccw around the current point/vector.
		Create a list for concavity test results.
			pv automatically fails.
			All other p are okay.
		For each p that has not failed:
			Cross the vectors from the current point to the previous okay point and the next okay point.
			Dot the result with pv.
			If the result is not positive, that point fails.
			Loop until for one lap after the last failure is discovered.
		For all the points that didn't fail:
			Create a triangle using pv, the point that has not failed, and the next (or previous) point.
	Return list of triangles.
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
		%disp(sprintf('ii = %2i/%2i', ii, len_p3))
		v = v3(ii,:);
		
		% get directional vectors
		% don't bother with the case v3(ii,:) == [0 0 0]
		% if that happens, something already went wrong
		u3 = v/norm(v);
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
			u1 = u1/norm(u1);
		elseif abs(v(2)) >= abs(v(1)) && abs(v(2)) >= abs(v(3))
			if abs(v(1)) >= abs(v(3))
				u1 = [v(2) -v(1) 0];
			else
				u1 = [0 v(3) -v(2)];
			end
			u1 = u1/norm(u1);
		else
			if abs(v(1)) >= abs(v(2))
				u1 = [v(3) 0 -v(1)];
			else
				u1 = [0 v(3) -v(2)];
			end
			u1 = u1/norm(u1);
		end
		u2 = cross(u3,u1);
		%disp(sprintf('u1 = [%f %f %f]', u1(1), u1(2), u1(3)))
		
		%disp('post unit vectors')
		
		% get angles of points around v
		theta = zeros(len_p3,1);
		phi   = zeros(len_p3,1);
		for jj = 1:len_p3
			theta(jj) = atan2( dot(u2,v3(jj,:)), dot(u1,v3(jj,:)) );
			phi(jj)   = dot(u3,v3(jj,:)) / norm(v3(jj,:));
		end
		
		% sort ccw
		[theta indices] = sort(theta, 'ascend');
		phi_s = phi(indices);
		ii_s = find(indices == ii, 1);
		v3_s = v3(indices,:);
		% inverse indices, 
		% so that they they can point to 
		% their counterparts on the unsorted lists
		[temp indexes] = sort(indices,'ascend');
		
		% test results, 1 means keep testing, 0 means failed
		t_s = ones(len_p3,1);
		% the origin of the fan need not apply
		t_s(ii_s) = 0;
	
		% cc keeps track of the position around the fan,
		% and the number of laps.
		% A more efficient method would be a single lap,
		% with backtracking whenever a concavity is discovered.
		% But for now, this should work.
		% This being that every time a concavity is found, 
		% the counter loops around again to ensure that 
		% the previous point is not also in the concavity.
		cc = 0;
		cc_max = len_p3+1;
		% jn is previous index (j-counter negative one)
		jn = mod(cc-1-1,len_p3)+1;
		while t_s(jn) == 0
			jn = jn-1;
			jn = mod(jn-1,len_p3)+1;
		end
		% jj is current index under test
		jj = mod(cc-1,len_p3)+1;
		while t_s(jj) == 0
			cc = cc+1;
			jj = mod(cc-1,len_p3)+1;
		end
		% jp is next index (j-counter positive one)
		% cc tracks this one because it's the leading edge
		cc = cc+1;
		jp = mod(cc-1,len_p3)+1;
		while t_s(jp) == 0
			cc = cc+1;
			jp = mod(cc-1,len_p3)+1;
		end
		while cc <= cc_max && sum(t_s) > 2
			% get points
			vn = v3_s(jn,:);
			vj = v3_s(jj,:); % My naming scheme went terribly wrong.
			vp = v3_s(jp,:);
		
			t = dot( cross(vn-v,vp-v), vj-v );
		
			%disp(sprintf('cc = %2i/%2i, jn = %2i, jj = %2i, jp = %2i, t = %+5f', cc, cc_max, jn, jj, jp, t))
		
			if t < 0
				t_s(jj) = 0;
				cc_max = cc + len_p3+1;
				jn = jn;
			else
				jn = jj;
			end
			jj = jp;
			% find nex valid point
			cc = cc+1;
			jp = mod(cc-1,len_p3)+1;
			while t_s(jp) == 0
				cc = cc+1;
				jp = mod(cc-1,len_p3)+1;
			end
		end
		
		%disp('post fan calc')
		
		%disp(t_s(indexes)')
		
		% create triangles
		cc = 0;
		cc_max = len_p3+1;
		jn = mod(cc-1-1,len_p3)+1;
		while t_s(jn) == 0
			jn = jn-1;
			jn = mod(jn-1,len_p3)+1;
		end		
		while cc <= cc_max
			jj = mod(cc-1,len_p3)+1;
			if t_s(jj) == 1
				% record triangle
				i2 = indices(jj);
				i3 = indices(jn);
				
				%{
				len_l3 = len_l3+1;
				l3(len_l3,:) = [ii i2 i3];
				%}
				if i2 < ii || i3 < ii
					% already accounted for
				elseif i2 < i3
					len_l3 = len_l3+1;
					l3(len_l3,:) = [ii i2 i3];
				else % if i2 < i3
					len_l3 = len_l3+1;
					l3(len_l3,:) = [ii i3 i2];
				end
				%}
				jn = jj;
			end
			cc = cc+1;
		end
		%disp('post fan creation')
	end
	
	% cull to proper size
	l3 = l3(1:len_l3,:);
	% sort triangles to make non-unique triangles stand out
	%l3 = transpose(sort(transpose(l3),'ascend'));
	% remove all non-unqiue triangles
	%l3 = unique(l3,'rows');
end
