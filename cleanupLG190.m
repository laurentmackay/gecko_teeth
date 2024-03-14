
function teethUpdated = cleanupLG190(teethOriginal,offset)

% add teeth (these were missed in the spreadsheet conversion somehow)
teethOriginal = [teethOriginal; 2 231 0; 3 11 0; 3 217 0; 4 116 0; 4 158  0; 4 200 0; 5 53 0; 16 39 0; 7 63 0; 8 11 0; 8 46 0;
    78 77 0; 78 161 0; 79 56 0; 80 74 0; 81 49 0];

teethOriginal=sortrows(teethOriginal,[1 2]);

% Remove teeth at (5,88), (59,84), (71,105), (72,126), (74,74), (79,74),
% (80,105) - these were incorrectly added in the conversion (e.g. 1 0 2 0 1)

ind = find(teethOriginal(:,1)==44 & teethOriginal(:,2)==4);
teethOriginal=teethOriginal([1:ind-1 ind+1:size(teethOriginal,1)],:);

ind = find(teethOriginal(:,1)==5 & teethOriginal(:,2)==88);
teethOriginal=teethOriginal([1:ind-1 ind+1:size(teethOriginal,1)],:);

ind = find(teethOriginal(:,1)==6 & teethOriginal(:,2)==21);
teethOriginal=teethOriginal([1:ind-1 ind+1:size(teethOriginal,1)],:);

ind = find(teethOriginal(:,1)==59 & teethOriginal(:,2)==84);
teethOriginal=teethOriginal([1:ind-1 ind+1:size(teethOriginal,1)],:);

ind = find(teethOriginal(:,1)==71 & teethOriginal(:,2)==105);
teethOriginal=teethOriginal([1:ind-1 ind+1:size(teethOriginal,1)],:);

ind = find(teethOriginal(:,1)==72 & teethOriginal(:,2)==126);
teethOriginal=teethOriginal([1:ind-1 ind+1:size(teethOriginal,1)],:);

ind = find(teethOriginal(:,1)==74 & teethOriginal(:,2)==74);
teethOriginal=teethOriginal([1:ind-1 ind+1:size(teethOriginal,1)],:);

ind = find(teethOriginal(:,1)==79 & teethOriginal(:,2)==74);
teethOriginal=teethOriginal([1:ind-1 ind+1:size(teethOriginal,1)],:);

ind = find(teethOriginal(:,1)==80 & teethOriginal(:,2)==105);
teethOriginal=teethOriginal([1:ind-1 ind+1:size(teethOriginal,1)],:);


% merge tooth numbers 65 and 66
teethOriginal(419:end,1)=teethOriginal(419:end,1)-1;

% move one eruption to test the plotting routine
% teeth(107,2) = teeth(107,2) + 13;

% Remove tooth numbers below 7 and above 76 where the bites seem to not
% always catch all teeth:
indKeep = find(abs(teethOriginal(:,1)-offset)<36.5);
teethOriginal = teethOriginal(indKeep,:);

teethUpdated = teethOriginal;

end
