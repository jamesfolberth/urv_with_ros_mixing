function [status] = add_lapack()
% A simple script to add lapack/ to the path or, if not found,
% download and install the package automatically.

status = 0;

if ~exist('lapack', 'dir') || ~exist('lapack/lapack.m', 'file')
   fprintf(['In order to run lapack.mex, you must first install the lapack interface\n'...
            'package from https://www.mathworks.com/matlabcentral/fileexchange/16777-lapack.\n'...
            'We can try to do this automatically for you.\n']);
   res = input('Proceed? (y/n) [y]:', 's'); 
   if isempty(res), res='y'; end
   
   if strcmp(lower(res), 'y')
      if ~exist('lapack', 'dir'), mkdir('lapack'); end

      fprintf('Downloading LAPACK interface package `lapack`.\n');
      URL = 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/16777/versions/9/download/zip';
      websave('lapack/lapack.zip', URL);
      unzip('lapack/lapack.zip', 'lapack/');
      addpath('lapack'); % for DGEQPX
      status = 1;
   else
      fprintf(['Please manually download the LAPACK interface package from\n'...
               'https://www.mathworks.com/matlabcentral/fileexchange/16777-lapack and\n'...
               'unzip and install the package in the newly created subdirectory\n'...
               '`lapack/`.\n']);
      return;
   end
else
   addpath('lapack'); % for DGEQPX
   status = 1;
end

end
