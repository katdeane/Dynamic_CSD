function [data] = chaninterp(data, method, chan_bad, pos_chan)

%channel interpolation first dimension of data
[nch, npt,nepoch,n] = size(data);


if ~exist('method','var') || isempty(method); method = 'linextra'; end;
if ~exist('chan_bad','var') || isempty(chan_bad); chan_bad = []; end;
if ~exist('pos_chan','var') || isempty(pos_chan); pos_chan = (1:nch); end;

nch_bad = length(chan_bad(:));
chan = (1:nch)';
chan_good = setdiff(chan, chan_bad);
%nch_good = length(chan_good(:)); %unused currently
nd = size(pos_chan,2);

if ~isempty(chan_bad)
    if isequal(lower(method),'linextra') 
        
        
            for ich=1:nch_bad
                

                chextra = [];
                chinter = [];
                if chan_bad(ich) == 1
                   chextra = chan_good(chan_good > chan_bad(ich));
                   if length(chextra(:)) >= 2; chextra = chextra(1:2); else warning('no extrapolation possible'); break; end;
                elseif chan_bad(ich) == max(chan(:));
                   chextra = chan_good(chan_good < chan_bad(ich));
                   if length(chextra(:)) >= 2; chextra = chextra(end-1:end); else warning('no extrapolation possible');break;  end;
                   chextra = chextra([2,1]);
                else
                   chinter1 = chan_good(chan_good < chan_bad(ich)); 
                   if length(chinter1(:)) >= 1; chinter1 = chinter1(end); else warning('no interpolation possible');break;  end;
                    chinter2 = chan_good(chan_good > chan_bad(ich)); 
                   if length(chinter2(:)) >= 1; chinter2 = chinter2(1); else warning('no interpolation possible');break;  end;
                   chinter = [chinter1, chinter2];
                end


                if ~isempty(chinter)
                   dchb = abs(chinter-chan_bad(ich));
                   data(chan_bad(ich),:,:) = ((1/dchb(1))*data(chinter(1),:,:) + (1/dchb(2))*data(chinter(2),:,:))/sum(1./dchb(:));
                end
                if ~isempty(chextra)
                   dchb  = abs(chextra-chan_bad(ich));
                   data(chan_bad(ich),:,:) = data(chextra(1),:,:) + (dchb(1)/diff(dchb)) * (data(chextra(1),:,:)-data(chextra(2),:,:));
                end
            end
    else 
        
        if nd ==1
            
            xbad = pos_chan(chan_bad,1);
            xgood = pos_chan(chan_good,1);
             for it = 1:npt
                    for ie = 1:nepoch
                        for in = 1:n

                            
                            [data(chan_bad,it,ie,in)] = interp1(xgood , data(chan_good,it,ie,in),xbad, method); % interpolate data               
                       
                                                               
                        end
                    end

             end
            %    yi = interp1(x,y,xi,method)
            %interp1
            %       Nearest neighbor interpolation (method = 'nearest'). This method sets the value of an interpolated point to the value of the nearest existing data point.
            %       Linear interpolation (method = 'linear'). This method fits a different linear function between each pair of existing data points, and returns the value of the relevant function at the points specified by xi. This is the default method for the interp1 function.
            %       Cubic spline interpolation (method = 'spline'). This method fits a different cubic function between each pair of existing data points, and uses the spline function to perform cubic spline interpolation at the data points.
            %       Cubic interpolation (method = 'pchip' or 'cubic'). These methods are identical. They use the pchip function to perform piecewise cubic Hermite interpolation within the vectors x and y. These methods preserve monotonicity and the shape of the data.
            % If any element of xi is outside the interval spanned by x, the specified interpolation method is used for extrapolation. Alternatively, yi = interp1(x,Y,xi,method,extrapval) replaces extrapolated values with extrapval. NaN is often used for extrapval.
            % All methods work with nonuniformly spaced data.
            % Speed, Memory, and Smoothness Considerations.   When choosing an interpolation method, keep in mind that some require more memory or longer computation time than others. However, you may need to trade off these resources to achieve the desired smoothness in the result:
            %       Nearest neighbor interpolation is the fastest method. However, it provides the worst results in terms of smoothness.
            %       Linear interpolation uses more memory than the nearest neighbor method, and requires slightly more execution time. Unlike nearest neighbor interpolation its results are continuous, but the slope changes at the vertex points.
            %       Cubic spline interpolation has the longest relative execution time, although it requires less memory than cubic interpolation. It produces the smoothest results of all the interpolation methods. You may obtain unexpected results, however, if your input data is nonuniform and some points are much closer together than others.
            %       Cubic interpolation requires more memory and execution time than either the nearest neighbor or linear methods. However, both the interpolated data and its derivative are continuous.

        elseif nd ==2
            xbad = pos_chan(chan_bad,1);
            ybad = pos_chan(chan_bad,2);
            xgood = pos_chan(chan_good,1);
            ygood = pos_chan(chan_good,2);
            for it = 1:npt
                    for ie = 1:nepoch
                        for in = 1:n
           %                 [Xi,Yi,data(chan_bad,it,ie,in)'] = griddata(ygood, xgood , data(chan_good,it,ie,in),...
           %                                                         ybad, xbad, method); % interpolate data              

                            [data(chan_bad,it,ie,in)] = interp2(xgood, ygood , data(chan_good,it,ie,in),...
                                                                    xbad, ybad, method); % interpolate data              
                                                                        % The function interp2 performs two-dimensional interpolation, an important operation for image processing and data visualization. Its most general form is
                % ZI = interp2(X,Y,Z,XI,YI,method)
                % Z is a rectangular array containing the values of a two-dimensional function, and X and Y are arrays of the same size containing the points for which the values in Z are given. XI and YI are matrices containing the points at which to interpolate the data. method is an optional string specifying an interpolation method.
                % There are three different interpolation methods for two-dimensional data:
                %       Nearest neighbor interpolation (method = 'nearest'). This method fits a piecewise constant surface through the data values. The value of an interpolated point is the value of the nearest point.
                %       Bilinear interpolation (method = 'linear'). This method fits a bilinear surface through existing data points. The value of an interpolated point is a combination of the values of the four closest points. This method is piecewise bilinear, and is faster and less memory-intensive than bicubic interpolation.
                %       Bicubic interpolation (method = 'cubic'). This method fits a bicubic surface through existing data points. The value of an interpolated point is a combination of the values of the sixteen closest points. This method is piecewise bicubic, and produces a much smoother surface than bilinear interpolation. This can be a key advantage for applications like image processing. Use bicubic interpolation when the interpolated data and its derivative must be continuous.
                % All of these methods require that X and Y be monotonic, that is, either always increasing or always decreasing from point to point. You should prepare these matrices using the meshgrid function, or else be sure that the "pattern" of the points emulates the output of meshgrid. In addition, each method automatically maps the input to an equally spaced domain before interpolating. If X and Y are already equally spaced, you can speed execution time by prepending an asterisk to the method string, for example, '*cubic'.



                        end
                    end

            end


        elseif nd ==3

            xbad = pos_chan(chan_bad,1);
            ybad = pos_chan(chan_bad,2);
            zbad = pos_chan(chan_bad,3);
            xgood = pos_chan(chan_good,1);
            ygood = pos_chan(chan_good,2);
            zgood = pos_chan(chan_good,3);

            for it = 1:npt
                    for ie = 1:nepoch
                        for in = 1:n
        %                    [data(chan_bad,it,ie,in)] = griddatan(pos_chan(chan_good,:), data(chan_good,it,ie,in),...
        %                                                            pos_chan(chan_bad,:), method); % interpolate data               

                         [data(chan_bad,it,ie,in)] = interp3(xgood,ygood,zgood, data(chan_good,it,ie,in),...
                                                                    xbad,ybad,zbad, method); % interpolate data              


                            %                                                         VI = interp3(X,Y,Z,V,XI,YI,ZI) interpolates to find VI, the values of the underlying three-dimensional function V at the points in arrays XI, YI and ZI. XI,YI, ZI must be arrays of the same size, or vectors. Vector arguments that are not the same size, and have mixed orientations (i.e. with both row and column vectors) are passed through meshgrid to create the Y1, Y2, Y3 arrays. Arrays X, Y, and Z specify the points at which the data V is given. Out of range values are returned as NaN.
                            % VI = interp3(V,XI,YI,ZI) assumes X=1:N, Y=1:M, Z=1:P where [M,N,P]=size(V).
                            % VI = interp3(V,ntimes) expands V by interleaving interpolates between every element, working recursively for ntimes iterations. The command interp3(V) is the same as interp3(V,1).
                            % VI = interp3(...,method) specifies alternative methods:
                            % 'nearest'
                            % Nearest neighbor interpolation
                            % 'linear'
                            % Linear interpolation (default)
                            % 'spline'
                            % Cubic spline interpolation
                            % 'cubic'
                            % Cubic interpolation, as long as data is uniformly-spaced. Otherwise, this method is the same as 'spline'

                        end
                    end

            end

        end
        
    end
end

    