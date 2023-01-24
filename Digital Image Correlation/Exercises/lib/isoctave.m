function a = isoctave
% return true when running octave
a = exist('OCTAVE_VERSION', 'builtin') ~= 0;
