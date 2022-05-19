warning('off','all');
% https://docs.google.com/presentation/d/1MmPZSztQQ9DxuBUmKQpQiRccXI-tjI7JR7H1n5OuJkU/edit?usp=sharing
% random_patch(sample, patchsize)
%     This method takes in a sample img and a 'patch' gets taken from
%     based on a random position. 
%     TLDR: Random subimage of sample.
%     
% quilt_random(sample, outsize, patchsize)
%     This method takes in a sample img creates mulitple random_patch's
%     and fills a outsize x outsize out img with them.
%     TLDR: Tiled random_patch's of size outsize x outsize.
% 
% quilt_simple(sample, outsizem patchsize, overlap, tol)
%     This method takes in sample then...
%        1. Places random_patch at top left.
%        2. Calculate ssd_patch between current patch position
%           and the left and above tile.
%        3. choose_sample based on ssd_patch where the error is
%           small based off tol.
%     TLDR:First tile random, rest try to share border similarities.
%
% ssd_patch(I,T,M)
%     This method performs template matching of the overlap and computes
%     the sum of squared differences.
%     TLDR:Calculate errors for overlapping region.
%     
% choose_sample(ssd, tol)
%     This method chooses a location where the ssd is low.
%     TLDR:Choose a small error.
%     
% quilt_cut(sample, outsizem patchsize, overlap, tol)
%     This method takes in sample then...
%        1. Places random_patch at top left.
%        2. Calculate ssd_patch between current patch position
%           and the left and above tile.
%        3. choose_sample based on ssd_patch where the error is
%           small based off tol.
%        4. Calculate shortest path between errors
%        5. Create a mask that cuts off border irregularities.
%        6. Add to result.
%     TLDR:First tile random, rest try to share border similarities and cut nonmatching borders.


code = imread('https://toppng.com/uploads/preview/code-symbols-text-programming-texture-11569700826bnlo2s3hxa.jpg');
random_code = quilt_random(code, 400, 40);
simple_code = quilt_simple(code, 400, 40, 5, 0.1);
cut_code = quilt_cut(code, 400, 40, 5, 0.1);
figure;montage({code,random_code,simple_code,cut_code});

cb = imread('https://i.ebayimg.com/images/g/IR0AAOSws-tbFceA/s-l400.jpg');
cut_cb = quilt_cut(cb, 500, 60, 5, 0.1);
figure;montage({cb,cut_cb});

fl = imread('https://images.squarespace-cdn.com/content/v1/561a8febe4b069644897cf4b/1484501601621-31B1QO72IY4BCG1QCVGE/DSC_0277.jpg');
cut_fl = quilt_cut(fl, 800,50, 5, 0.1);
figure;montage({dl,cut_fl});

wp = imread('https://garden.spoonflower.com/c/10064285/p/f/m/xMD6e3QC1Xl3RRMt-bCZzzvuTTo82jzkQ56O-YgfqNw2gCLpEBWW1bk/paint%20splatter%20repeat%201.jpg');
cut_wp = quilt_cut(wp, 400,40, 8, 0.1);
figure;montage({wp,cut_wp});


gr = imread('https://drive.google.com/uc?id=1jc9Fzm4WymRYB6IMpI2Lxd9x1dmPutfL', 'jpg');
gr = gr(1:1512,1:1512,:);
cut_gr = quilt_cut(gr, 750, 50, 8, 0.1);
figure;montage({gr,cut_gr});

function bool = showborders()
    % 1 = show the cut borders
    % 0 = don't show the cut borders
    bool = 1;
end

function imout = quilt_random(sample, outsize, patchsize)
    for i = 1:patchsize:outsize
        for j = 1:patchsize:outsize
            patch = random_patch(sample, patchsize);  
            imout(i:i+patchsize-1,j:j+patchsize-1,:) = patch;
        end
    end
end

function imout = quilt_simple(sample, outsize, patchsize, overlap, tol)
    left = zeros(patchsize,patchsize);
    left(:,1:overlap) = 1;

    top = zeros(patchsize,patchsize);
    top(1:overlap,:) = 1;

    topleft = top + left;
    topleft(topleft>1)=1;

    imout = zeros(outsize, outsize, 3);
    imout = uint8(imout);

    numBlocks = (outsize-overlap) / (patchsize-overlap);

    for i = 1:numBlocks
        for j = 1:numBlocks
            if i == 1 && j == 1
                patch = random_patch(sample, patchsize);
                imout(1:patchsize, 1:patchsize,:) = patch;
                continue
            elseif i == 1
                M = left;
                I = sample;
                y = (i-1)*(patchsize-overlap) + 1;
                x = (j-1)*(patchsize-overlap) + 1;
                T = imout(y:y+patchsize-1, x:x+patchsize-1, :);
                ssd = ssd_patch(I, T, M);  

                [r, c] = choose_sample(ssd, tol);
                patch = sample(r:r+patchsize-1, c:c+patchsize-1, :);
            elseif j == 1
                M = top;
                I = sample;
                y = (i-1)*(patchsize-overlap) + 1;
                x = (j-1)*(patchsize-overlap) + 1;
                T = imout(y:y+patchsize-1, x:x+patchsize-1, :);
                ssd = ssd_patch(I, T, M);  

                [r, c] = choose_sample(ssd, tol);
                patch = sample(r:r+patchsize-1, c:c+patchsize-1, :);
            else
                M = topleft;
                I = sample;
                y = (i-1)*(patchsize-overlap) + 1;
                x = (j-1)*(patchsize-overlap) + 1;
                T = imout(y:y+patchsize-1, x:x+patchsize-1, :);
                ssd = ssd_patch(I, T, M);  

                [r, c] = choose_sample(ssd, tol);
                patch = sample(r:r+patchsize-1, c:c+patchsize-1, :);
            end
            startPos = [(i-1)*(patchsize-overlap)+1, (j-1)*(patchsize-overlap)+1];
            endPos =   [(i-1)*(patchsize-overlap)+patchsize, (j-1)*(patchsize-overlap)+patchsize];
            imout(startPos(1):endPos(1),startPos(2):endPos(2),:) = patch;
        end
    end
end

function imout = quilt_cut(sample, outsize, patchsize, overlap, tol)
    left = zeros(patchsize,patchsize);
    left(:,1:overlap) = 1;

    top = zeros(patchsize,patchsize);
    top(1:overlap,:) = 1;

    topleft = top + left;
    topleft(topleft>1)=1;
    figure;montage({left, top, topleft});

    imout = zeros(outsize, outsize, 3);
    imout = uint8(imout);

    numBlocks = (outsize-overlap) / (patchsize-overlap);

    for i = 1:numBlocks
        for j = 1:numBlocks
            if i == 1 && j == 1
                patch = random_patch(sample, patchsize);
                imout(1:patchsize, 1:patchsize,:) = patch;
                continue
            elseif i == 1
                M = left;
                I = sample;
                y = (i-1)*(patchsize-overlap) + 1;
                x = (j-1)*(patchsize-overlap) + 1;
                T = imout(y:y+patchsize-1, x:x+patchsize-1, :);
                ssd = ssd_patch(I, T, M);  

                [r, c] = choose_sample(ssd, tol);
                patch = sample(r:r+patchsize-1, c:c+patchsize-1, :);

                patch1 = imout(1:patchsize, x:x+overlap-1,:);
                patch2 = patch(1:patchsize, 1:overlap, :);

                err = sum((im2double(patch1) - im2double(patch2)).^2, 3);
                mask = transpose(cut(transpose(err)));
                mask = [mask ones(patchsize, patchsize-overlap)];
            elseif j == 1
                M = top;
                I = sample;
                y = (i-1)*(patchsize-overlap) + 1;
                x = (j-1)*(patchsize-overlap) + 1;
                T = imout(y:y+patchsize-1, x:x+patchsize-1, :);
                ssd = ssd_patch(I, T, M);  

                [r, c] = choose_sample(ssd, tol);
                patch = sample(r:r+patchsize-1, c:c+patchsize-1, :);

                patch1 = imout(y:y+overlap-1, 1:patchsize,:);
                patch2 = patch(1:overlap, 1:patchsize, :);

                err = sum((im2double(patch1) - im2double(patch2)).^2, 3);
                mask = cut(err);
                mask = [mask;ones(patchsize-overlap, patchsize)];
            else
                M = topleft;
                I = sample;
                y = (i-1)*(patchsize-overlap) + 1;
                x = (j-1)*(patchsize-overlap) + 1;
                T = imout(y:y+patchsize-1, x:x+patchsize-1, :);
                ssd = ssd_patch(I, T, M);  

                [r, c] = choose_sample(ssd, tol);
                patch = sample(r:r+patchsize-1, c:c+patchsize-1, :);

                patch1 = imout(y:y+patchsize-1, x:x+overlap-1,:);
                patch2 = patch(1:patchsize, 1:overlap, :);

                err = sum((im2double(patch1) - im2double(patch2)).^2, 3);
                mask = transpose(cut(transpose(err)));
                mask1 = [mask ones(patchsize, patchsize-overlap)];



                patch1 = imout(y:y+overlap-1, x:x+patchsize-1,:);
                patch2 = patch(1:overlap, 1:patchsize, :);

                err = sum((im2double(patch1) - im2double(patch2)).^2, 3);
                mask = cut(err);
                mask2 = [mask;ones(patchsize-overlap, patchsize)];

                mask = mask1 & mask2;
            end
            
            mask = uint8(mask);
            startPos = [(i-1)*(patchsize-overlap)+1, (j-1)*(patchsize-overlap)+1];
            endPos =   [(i-1)*(patchsize-overlap)+patchsize, (j-1)*(patchsize-overlap)+patchsize];
            patch_keep = imout(startPos(1):endPos(1),startPos(2):endPos(2),:) .* uint8(~mask);
            patch_cut = uint8(patch) .* mask;
            if showborders()
                boundary = boundarymask(mask);
                patch_cut(boundary==1) = 0; 
            end
            if showborders() && i==1 && j==2
                figure;montage({patch_cut,boundary});
            end
            
            
            imout(startPos(1):endPos(1),startPos(2):endPos(2),:) = patch_cut + patch_keep;
        end
    end
end

function ssd = ssd_patch(I, T, M)
    I = im2double(I);
    T = im2double(T);
    M = im2double(M);
    A = I;
    A(:,:,1) = 2*imfilter(A(:,:,1), M.*T(:,:,1));
    A(:,:,2) = 2*imfilter(A(:,:,2), M.*T(:,:,2));
    A(:,:,3) = 2*imfilter(A(:,:,3), M.*T(:,:,3));
    ssd = imfilter(I.^2, M) -A + sum(sum((M.*T).^2));
    ssd = sum(ssd,3);
    Isize = size(ssd);
    patchsize = size(T,1);
    ssd = ssd((1+patchsize)/2:Isize(1)-(1+patchsize)/2+1, (1+patchsize)/2:Isize(2)-(1+patchsize)/2+1, :);
end

function [y, x] = choose_sample(ssd, tol)
minc = min(min(ssd));
minc=max(minc,.0001);
[r, c] = find(ssd<=minc*(1+tol));
x = randi([1 size(r,1)],1,1);
y = r(uint8(x));
x = c(uint8(x));
end

function mask = cut(errpatch)
% mask = cut(errpatch)
%
% Computes the minimum cut path from the left to right side of the patch
% 
% Input:
%   errpatch: cost of cutting through each pixel
% Output:
%   mask: a 0-1 mask that indicates which pixels should be on either side
%   of the cut


% create padding on top and bottom with very large cost 
errpatch = [1E10*ones(1, size(errpatch,2)) ; errpatch ; 1E10*ones(1, size(errpatch,2))];
[h, w] = size(errpatch);

path = zeros([h w]);

cost = zeros([h w]);
cost(:, 1) = errpatch(:, 1);
cost([1 end], :) = errpatch([1 end], :);
for x = 2:w  % for each column, compute the cheapest connected path to the left
  % cost of path for each row from left upper/same/lower pixe
  tmp = [cost(1:h-2, x-1) cost(2:h-1, x-1) cost(3:h, x-1)]; 
  [mv, mi] = min(tmp, [], 2);  % mi corresponds to upper/same/lower for each row
  path(2:h-1, x) = (2:(h-1))'+mi-2; % save the next step of the path
  cost(2:h-1, x) = cost(path(2:h-1, x), x-1)+errpatch(2:h-1, x);  % update the minimum costs for each 
end
path = path(2:end-1, :)-1; % remove padding
cost = cost(2:end-1, :);

% create the mask based on the best path
mask = zeros(size(path));
bestpath = zeros(1, size(path, 2));
[mv, bestpath(end)] = min(cost(:, end)); 
mask(1:bestpath(end), end) = 1;

for x = numel(bestpath):-1:2
  bestpath(x-1) = path(bestpath(x), x);
  mask(1:bestpath(x-1), x-1) = 1;
end
mask = ~mask;
end

function im = random_patch(sample, patchsize)
    [y,x,n] = size(sample);
    yStart = randi([1 y-patchsize],1,1);
    xStart = randi([1 x-patchsize],1,1);
    im = sample(yStart:yStart+patchsize-1,xStart:xStart+patchsize-1,:);
end