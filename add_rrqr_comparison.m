function [status] = add_rrqr_comparison()
% A simple script to add rrqr_comparision/ to the path or, if not found,
% download and install the package automatically.

status = 0;

%NOTE: In order to run DGEQPX, you must first install the rank revealing codes
%      package from http://www.math.sjsu.edu/~foster/rankrevealingcode.html.
if ~exist('rrqr_comparison', 'dir') || ~exist('rrqr_comparison/qrxp.m', 'file')
   fprintf(['In order to run DGEQPX, you must first install the rank revealing\n'...
            'codes package from http://www.math.sjsu.edu/~foster/rankrevealingcode.html.\n'...
            'We can try to do this automatically for you.\n']);
            %'Please follow the link and install the appropriate package in the\n'...
            %'newly created subdirectory `rrqr_comparison` and re-run qlp_test\n']);
   res = input('Proceed? (y/n) [y]:', 's'); 
   if isempty(res), res='y'; end
   
   if strcmp(lower(res), 'y')
      if ~exist('rrqr_comparison', 'dir'), mkdir('rrqr_comparison'); end

      res = input(['The rank revealing codes package contains pre-compiled mex\n'...
                   'library files for\n'...
                   ' [0] - for Linux computers with 64 bit Matlab and an Intel x86-64 processor\n'...
                   ' [1] - for Windows computers with 64 bit Matlab 7.6 or higher\n'...
                   ' [2] - for Windows computers with 32 bit Matlab 7.6 or higher\n'...
                   ' [-1] - skip automatic installation\n'...
                   'Please select an available architecture [-1]:'], 's');
      res = str2num(res);
      if isempty(res), res=-1; end

      if res == -1
         fprintf(['Please manually download the appropriate package from\n'...
                  'http://www.math.sjsu.edu/~foster/rankrevealingcode.html and\n'...
                  'unzip and install the package in the newly created subdirectory\n'...
                  '`rrqr_comparision/`.\n']);
         return
      elseif res == 0
         URL = 'http://www.math.sjsu.edu/~foster/rank/rank_matlab_linux_64bit.zip';
      elseif res == 1
         URL = 'http://www.math.sjsu.edu/~foster/rank/rank_matlab_7p6_7_8_64bit.zip';
      elseif res == 2
         URL = 'http://www.math.sjsu.edu/~foster/rank/rank_matlab_7p6_7_8_32bit.zip';
      else
         return
      end
      
      fprintf('Downloading rank revealing codes package.\n');
      websave('rrqr_comparison/rank.zip', URL);
      unzip('rrqr_comparison/rank.zip', 'rrqr_comparison/');
      addpath('rrqr_comparison'); % for DGEQPX
      status = 1;
   else
      fprintf(['Please manually download the appropriate package from\n'...
               'http://www.math.sjsu.edu/~foster/rankrevealingcode.html and\n'...
               'unzip and install the package in the newly created subdirectory\n'...
               '`rrqr_comparision/`.\n']);
      return;
   end
else
   addpath('rrqr_comparison'); % for DGEQPX
   status = 1;
end

end
