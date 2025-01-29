function [C, weights] = Crealdata(directory)

% Strain measurements
calib = load(directory + "strain0_sig.mat").results;
% valid = load(directory + "strain90_sig.mat").results;
idx = [10 17 24 28]; % AA positions

% Curvatures
kcalib = [0.5; 1.6; 2; 2.5; 3.2; 4];
% kvalid = [0.25; 0.8; 1; 1.25; 3.125];
angcalib = 0:20:340;
% angvalid = 10:20:350;
% kx = [kcalib*cosd(angcalib); kvalid*cosd(angvalid)]';
% ky = [kcalib*sind(angcalib); kvalid*sind(angvalid)]';
kx = (kcalib*cosd(angcalib))';
ky = (kcalib*sind(angcalib))';
k = [kx(:), ky(:)];

% Strain
[sig1, sig2, sig3] = ExtractSig(calib, kcalib, idx);
% [sig1v, sig2v, sig3v] = ExtractSig(valid, kvalid, idx);
% sig1 = [sig1c; sig1v];
% sig2 = [sig2c; sig2v];
% sig3 = [sig3c; sig3v];

% Temperature compensation
[AA1, AA2, AA3, AA4] = AATcomp(sig1, sig2, sig3);

% Linear mapping
CAA1 = pinv(AA1)*k;
CAA2 = pinv(AA2)*k;
CAA3 = pinv(AA3)*k;
CAA4 = pinv(AA4)*k;
C = cat(3, CAA1, CAA2, CAA3, CAA4);

% Fitting validation
KAA1 = AA1*CAA1;
KAA2 = AA2*CAA2;
KAA3 = AA3*CAA3;
KAA4 = AA4*CAA4;

% Errors
MSE1 = mean((k - KAA1).^2, "all");
MSE2 = mean((k - KAA2).^2, "all");
MSE3 = mean((k - KAA3).^2, "all");
MSE4 = mean((k - KAA4).^2, "all");
MSE = [MSE1 MSE2 MSE3 MSE4];

% Reliability weights for each sensing area
weights = 1./MSE;
weights = weights/sum(weights); % Normalize so sum = 1

% Weighted shape
KAA1w = weights(1)*KAA1;
KAA2w = weights(2)*KAA2;
KAA3w = weights(3)*KAA3;
KAA4w = weights(4)*KAA4;
k_weight = KAA1w + KAA2w + KAA3w + KAA4w;

% Error table
err1 = [mean(abs(k-KAA1), "all"), std(abs(k-KAA1), 0, "all")];
err2 = [mean(abs(k-KAA2), "all"), std(abs(k-KAA2), 0, "all")];
err3 = [mean(abs(k-KAA3), "all"), std(abs(k-KAA3), 0, "all")];
err4 = [mean(abs(k-KAA4), "all"), std(abs(k-KAA4), 0, "all")];
err = [err1; err2; err3; err4];

% Weighted error
errw = [mean(abs(k-k_weight), "all"), std(abs(k-k_weight), 0, "all")];

%%%%%%%%
function [sig1, sig2, sig3] = ExtractSig(data, k, idx)

    sig1 = [];
    sig2 = [];
    sig3 = [];

    for i = 1:length(k)
    
        field = fieldnames(data.fiber1{1});
    
        data1 = data.fiber1{1}.(field{i});
        sig1 = [sig1; data1(:,idx)];
        data2 = data.fiber2{1}.(field{i});
        sig2 = [sig2; data2(:,idx)];
        data3 = data.fiber3{1}.(field{i});
        sig3 = [sig3; data3(:,idx)];
    
    end

end

function [AA1, AA2, AA3, AA4] = AATcomp(sig1, sig2, sig3)

    AA1 = [sig1(:,1) sig2(:,1) sig3(:,1)];
    AA2 = [sig1(:,2) sig2(:,2) sig3(:,2)];
    AA3 = [sig1(:,3) sig2(:,3) sig3(:,3)];
    AA4 = [sig1(:,4) sig2(:,4) sig3(:,4)];
    
    AA1 = AA1 - mean(AA1,2);
    AA2 = AA2 - mean(AA2,2);
    AA3 = AA3 - mean(AA3,2);
    AA4 = AA4 - mean(AA4,2);

end
%%%%%%%%

end % function Crealdata