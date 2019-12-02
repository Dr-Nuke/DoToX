            im=im2double(fitsread(strcat(rpath,rfile))');

            % crop away area of bad DC image
            % border pixel calculation

            im=im(81:end,201:2250);
size(im)
size(im(6:end-5,6:end-5))