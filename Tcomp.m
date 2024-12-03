function [wvl_shift] = Tcomp(wvl_shift)

Tshift = mean(wvl_shift, 1);
wvl_shift = wvl_shift - Tshift;

end % function Tcomp