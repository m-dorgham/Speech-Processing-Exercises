function [DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogramA, mLogMelSpectrogramB)
	for i=1:columns(mLogMelSpectrogramA)
		for j=1:columns(mLogMelSpectrogramB)
			mDistance(i,j) = sqrt(sum((mLogMelSpectrogramA(:,i)-mLogMelSpectrogramB(:,j)).^2));
		end
	end
	%initialize mAccDistance with the cumulative sum of mDistance
	mAccDistance = cumsum(mDistance); %accumulates over columns
	mAccDistance(1,:) = cumsum(mDistance(1,:));
	mPrevious = zeros(rows(mAccDistance), columns(mAccDistance));
	%initialize mPrevious first row elements to be of horizontal direction
	mPrevious(1,2:end) = 1;
	%initialize mPrevious first column elements to be of vertical direction
	mPrevious(2:end,1) = 3;

	for n=2:rows(mAccDistance)
		for m=2:columns(mAccDistance)
			% compute accumulated distance for horizontal , diagonal and vertical step
			Horizontal = mAccDistance(n, m-1) + mDistance(n, m);
			Diagonal = mAccDistance(n-1, m-1) + 2*mDistance(n, m);
			Vertical = mAccDistance(n-1, m) + mDistance(n, m);
			% compute the minimum cost and remember the direction of the optimal path.
			[mAccDistance(n, m), mPrevious(n, m)] = min([Horizontal, Diagonal, Vertical]);
		end
	end

	%normalize to obtain the DTW cost
	DTWCost = mAccDistance(end,end) / ( rows(mAccDistance) + columns(mAccDistance) );

	%construct the path
	n = rows(mAccDistance);
	m = columns(mAccDistance);
	prev = mPrevious(n, m);
	vPath(1) = prev;
	idx = 2;
	while prev > 0
		if prev == 1
			vPath(idx) = mPrevious(n, m-1);
			m = m-1;
		elseif prev == 2
			vPath(idx) = mPrevious(n-1, m-1);
			n = n-1;
			m = m-1;
		elseif prev == 3
			vPath(idx) = mPrevious(n-1, m);
			n = n-1;
		end
		prev = vPath(idx);
		++idx;
	end
endfunction


