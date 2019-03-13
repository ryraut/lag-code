function FORMAT = create_blocks(format,min_length,tr)
	FORMAT = {};
	curr_block=[];
	j = 1;
	for i = 1:numel(format)
		if format(i)
			curr_block = [curr_block i]; %#ok<*AGROW>
		else
			if numel(curr_block)*tr >= min_length
				FORMAT{j} = curr_block;
				j = j + 1; 
			end
			curr_block=[];
		end
	end
	if numel(curr_block)*tr >= min_length
		FORMAT{j} = curr_block;
	end

end
