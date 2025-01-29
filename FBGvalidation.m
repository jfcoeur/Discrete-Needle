function [s_FBG, kappa_FBG] = FBGvalidation(directory, C)

% FBG positions
s_FBG = load(directory + "s_FBG.mat").s_FBG;

% Validation strain data
valid = load(directory + "strain90_sig.mat").results;
idx = [10 17 24 28]; % AA positions
[sig1, sig2, sig3] = ExtractSig(valid, idx);
sig = cat(3, sig1, sig2, sig3);

% Temperature compensation
sig = Tcomp(sig);

% FBG experimental curvature
for i = 1:size(sig,2)

    AA = reshape(sig(:,i,:),[size(sig,1),3]);
    kappa_FBG(:,:,i) = AA*C(:,:,i);

end

kappa_FBG(:,3,:) = zeros(size(kappa_FBG(:,1,:)));


%%%%%%%%
function [sig1, sig2, sig3] = ExtractSig(data, idx)

    sig1 = [];
    sig2 = [];
    sig3 = [];

    for j = 1:5
    
        field = fieldnames(data.fiber1{1});
    
        data1 = data.fiber1{1}.(field{j});
        sig1 = [sig1; data1(:,idx)];
        data2 = data.fiber2{1}.(field{j});
        sig2 = [sig2; data2(:,idx)];
        data3 = data.fiber3{1}.(field{j});
        sig3 = [sig3; data3(:,idx)];
    
    end

end
%%%%%%%%

end % function FBGcurvature